#! /usr/bin/env ruby


########################################################
DIR = File.dirname(__FILE__)
$: << File.join(DIR, 'lib')
$: << File.expand_path("~/project/Rhizobiales/scripts/")


########################################################
require 'getoptlong'
require 'parallel'

require 'Dir'
require 'util'
require 'chang_yong'
require 'processbar'
require 'extractCOG'
require 'seqIO'
require 'kegg'


########################################################
$CDD_TBL = File.expand_path("~/resource/db/cdd/cddid_all.tbl")
$CDD_BITSCORE_FILE = File.expand_path("~/resource/db/cdd/bitscore_specific_3.16.txt")
$CDD_SUPER_RELA_FILE = File.expand_path("~/resource/db/cdd/family_superfamily_links")

$KEGG_DIR = File.expand_path("~/resource/db/kegg/latest")
#$KEGG_GENE_2_ORTHO = File.join($KEGG_DIR, 'map.from_gene_to_orthology')
$KEGG_GENE_2_ORTHOs = Array.new
$KEGG_GENE_2_ORTHOs << File.join($KEGG_DIR, 'kegg_bacteria.gene2ko')

#$MODEL_ORGNS = %w[Bradyrhizobium_elkanii_USDA_76 Mesorhizobium_loti_MAFF303099 Rhizobium_etli_CFN_42 Sinorhizobium_meliloti_1021 Bartonella_henselae_U4 Brucella_abortus_2308]


########################################################
indirs = Array.new
outdir = nil
species_list_file = nil
seq_indirs = Array.new
seq_suffix = 'protein'
is_output_seq = false
seq_outdir = nil
fam_include_list_files = Array.new
type = nil
evalue_cutoff = 1e-3
cdd_type = nil
is_bitscore = false
cpu = 1
locus_outdir = nil
is_super_fam = false
is_force = false
is_tolerate = false
is_best = true


cdd_objs = Hash.new
gene2ortho = Hash.new

taxa_included = Hash.new
infiles = Array.new
seq_infiles = Array.new
seq_objs = Hash.new
fams_included = Array.new
pfams = Hash.new
fam2super = Hash.new


########################################################
class String
  def famModify
    str = self
    if self =~ /^(PF[^.]+).+/
      str = sub($&, $1).downcase.sub('pf', 'pfam')
    elsif self =~ /^ko:/
      str = sub(/^ko:/, '')
    end
    return(str)
  end
end


class Pfam
  attr_accessor :id, :name, :taxaInfo, :evalue, :start, :stop
  def initialize(pfam_id)
    @id = pfam_id
    @taxaInfo = Hash.new
  end
  def isOverlap?(given_start, given_stop)
    if given_start > @stop or given_stop < @start
      return(false)
    else
      return(true)
    end
  end
end


class Gene
  attr_accessor :name, :pfams
  def initialize(name)
    @name = name
    @pfams = Array.new
  end
end


########################################################
def getTypeBasedOnCddType(type, cdd_type)
  if cdd_type.nil?
    STDERR.puts "cdd_type is not given! Exiting ......"
    exit 1
  end
  if not type.nil?
    if type !~ /^pfam|hmmsearch|blast|kegg|kofam$/
      STDERR.puts "type #{type} wrong! Exiting ......"
      exit 1
    end
  else
    case cdd_type
      when /CDD|CL|SMART|COG|KOG|TIGRFAM|PRK/
        type = 'blast'
      when /hmmsearch|pfam/
        type = 'pfam'
      when /kegg/
        type = 'kegg'
      when /kofam/
        type = 'kofam'
    end
  end
  return(type)
end


def getSuperFamRela(infile)
  fam2super = Hash.new
  in_fh = File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    fam_name, fam_id, cl_name, cl_id = line_arr
    fam2super[fam_id] = cl_id
  end
  in_fh.close
  return(fam2super)
end


def filterCddType(cdd_objs, cdd_type)
  cdd_types = [:CDD, :CL, :SMART, :COG, :KOG, :TIGRFAM, :PRK]
  cdd_type_prefix_info = {
    :CDD => 'cdd',
    :CL => 'cl',
    :COG => 'COG',
    :KOG => 'KOG',
    :TIGRFAM => 'TIGR',
    :SMART => 'smart',
    :PRK => 'PRK',
  }

  #unless cdd_type_prefix_info.include?(cdd_type)
  #  puts "cdd_type #{cdd_type} wrong! Exiting ......"; exit 1
  #end

  new_cdd_objs = Hash.new

  cdd_objs.each_pair do |cdd, cdd_obj|
    if cdd_obj.accn =~ /^([A-Z]+)\d+/i
      if cdd_type_prefix_info[cdd_type] == 'cl'
        new_cdd_objs[cdd] = cdd_obj
      elsif cdd_type_prefix_info[cdd_type] == $1
        new_cdd_objs[cdd] = cdd_obj
      end
    end
  end

  return(new_cdd_objs)
end


def readSearchResult(type, infile, evalue_cutoff, cdd_objs, gene2ortho, fam2super, is_bitscore, is_best)
  genes = Hash.new
  pfams = Hash.new
  in_fh = File.open(infile, 'r')
  taxon = getCorename(infile)

  in_fh.each_line do |line|
    line.chomp!
    next if line =~ /^#/
    line_arr = line.split(/\s+/)
    #gene_name, pfam_name, pfam_id, qstart, qstop, evalue = 1.upto(6).map{|i|nil}

    case type
      when /pfam|hmmsearch/
        #Xanthobacter_autotrophicus_Py2|XAUT_RS11885 -            212 1-cysPrx_C           PF10417.7     40   2.5e-15   54.2   0.1   1   1   1.8e-18   4.5e-15   53.5   0.1     1    40   157   196   157   196 0.97 [NC_009720:c(2637921..2638559)] [peroxidase]
        gene_name, pfam_name, pfam_id = line_arr.values_at(0,3,4)
        qstart, qstop = line_arr.values_at(17,18).map{|i|i.to_i}
        evalue = line_arr[12].to_f
      when 'blast'
        #gene_name	gnl|CDD|223104	35.28	360	208	7	76	425	80	424	8e-60	  206
        gene_name = line_arr[0]
        cdd = line_arr[1].split(/[|:]/)[-1]
        cdd = fam2super[cdd] if not fam2super.empty?
        next if not cdd_objs.include?(cdd)
        cdd_obj = cdd_objs[cdd]
        bit_score = line_arr[-1].to_f
        next if bit_score < cdd_obj.bit_score if is_bitscore
        pfam_id = cdd_obj.accn
        evalue = line_arr[10].to_f
        qstart, qstop = line_arr[6,2].map{|i|i.to_i}
      when 'kegg'
        gene_name = line_arr[0]
        subject = line_arr[1]
        next if not gene2ortho.include?(subject)
        pfam_id = gene2ortho[subject]
        evalue = line_arr[10].to_f
        qstart, qstop = line_arr[6,2].map{|i|i.to_i}
      when 'kofam'
        #Azorhizobium_caulinodans_ORS_571|AZC_RS05320 K02587  501.43  857.7  6.7e-261 nitrogenase molybdenum-cofactor synthesis protein NifE
        gene_name, pfam_id = line_arr[1, 2] # note that the first item is ""
        pfam_name = pfam_id
        evalue = line_arr[5].to_f
        qstart, qstop = 1, 2
    end

    pfam = Pfam.new(pfam_id)
    pfam.name = pfam_name
    pfam.evalue = evalue
    if not qstart.nil? and not qstop.nil?
      pfam.start = qstart
      pfam.stop = qstop
    end

    if evalue <= evalue_cutoff
      if not genes.include?(gene_name)
        genes[gene_name] = Gene.new(gene_name)
      end
      gene = genes[gene_name]
      gene.pfams << pfam
    end
  end

  genes.each_pair do |gene_name, gene|
    good_pfams = is_best ? select_good_pfams(gene.pfams) : gene.pfams
    good_pfams.each do |pfam|
      pfam = ! pfams.has_key?(pfam.id) ? Pfam.new(pfam.id) : pfams[pfam.id]
      if not pfam.taxaInfo.include?(taxon.to_sym)
        pfam.taxaInfo[taxon.to_sym] = Array.new
      end
      pfam.taxaInfo[taxon.to_sym] << gene_name
      pfams[pfam.id] = pfam
    end
  end

  in_fh.close

  return(pfams)
end


def select_good_pfams(pfams)
  good_pfams = Array.new
  pfams.sort_by{|pfam|pfam.evalue}.each do |pfam|
    if good_pfams.empty?
      good_pfams << pfam
      break if pfam.start.nil?
    else
      start, stop = pfam.start, pfam.stop
      if good_pfams.map{|i|i.isOverlap?(start, stop)}.all?{|i|i == false}
        good_pfams << pfam
      end
    end
  end
  #p good_pfams.map{|i|[i.start, i.stop]} if good_pfams.size >= 2
  return(good_pfams)
end


def output_locus(pfams, locus_outdir)
  taxon2outfh = Hash.new
  pfams.each_pair do |pfam_id, pfam|
    pfam.taxaInfo.each_pair do |taxon, genes|
      outfile = File.join(locus_outdir, taxon.to_s)
      out_fh = taxon2outfh.include?(taxon) ? taxon2outfh[taxon] : File.open(File.join(locus_outdir, taxon.to_s), 'a')
      out_fh.puts [pfam_id, genes].flatten.join("\t")
    end
  end
  taxon2outfh.each_value do |out_fh|
    out_fh.close
  end
  puts "Locus output done!" or exit 0
end


def runInParallel(type, infiles, cpu, evalue_cutoff, cdd_objs, gene2ortho, seq_objs, fams_included, fam2super, seq_outdir, locus_outdir, taxa_included, is_bitscore, is_best)
  pfams = Hash.new

  results = Parallel.map(infiles, in_threads: cpu) do |infile|
    pfams = readSearchResult(type, infile, evalue_cutoff, cdd_objs, gene2ortho, fam2super, is_bitscore, is_best)
    pfams
  end

  results.each do |i|
    i.each_pair do |pfam_id, pfam|
      pfams[pfam_id] = pfam if not pfams.has_key?(pfam_id)
      pfams[pfam_id].taxaInfo.merge!(pfam.taxaInfo)
    end
  end

  output_locus(pfams, locus_outdir) unless locus_outdir.nil?
  output_pfams(pfams, fams_included, seq_objs, seq_outdir, type, taxa_included) if not seq_outdir.nil?
  return(pfams)
end


def output_pfams(pfams, fams_included, seq_objs, seq_outdir, type, taxa_included)
  mkdir_with_force(File.join(seq_outdir,'seq'), true)
  info_outfile = File.join(seq_outdir, 'seq_info.tbl')
  info_out_fh = File.open(info_outfile, 'w')
  pfams, new2oldPfam = create_new_pfams(pfams, type)

  all_taxa = fams_included.select{|pfam_name|pfams.include?(pfam_name)}.map{|pfam_name|pfams[pfam_name].taxaInfo.keys}.flatten.uniq

  fams_included.each do |pfam_name|
    next if not pfams.include?(pfam_name)
    info_output_items = Array.new
    old_pfam_name = new2oldPfam[pfam_name]
    seq_outfile = File.join(seq_outdir, 'seq', old_pfam_name + '.fas')
    seq_out_fh = File.open(seq_outfile, 'w')

    info_output_items << pfam_name
    pfam = pfams[pfam_name]

    taxa = taxa_included.empty? ? all_taxa : taxa_included.keys
    taxa.select{|taxon|seq_objs.include?(taxon.to_s)}.each do |taxon|
      if pfam.taxaInfo.include?(taxon.to_sym)
        ele = [taxon, pfam.taxaInfo[taxon.to_sym].join(',')].join(':')
        pfam.taxaInfo[taxon.to_sym].each do |gene|
          seq_out_fh.puts '>' + gene
          seq_out_fh.puts seq_objs[taxon.to_s][gene].seq
        end
      else
        ele = ''
      end
      info_output_items << ele
    end

    seq_out_fh.close
    info_out_fh.puts info_output_items.join("\t")
  end
  info_out_fh.close
end


def outputResult(pfams, type, outdir)
  count = 0
  STDOUT.puts "Outputting results ......"
  pfams.each_pair do |pfam_id, pfam|
    if type == 'kegg' 
      outfile = File.join(outdir, pfam_id.sub(/^(ko:)?/, ''))
    else
      outfile = File.join(outdir, pfam_id)
    end
    out_fh = File.open(outfile, 'w')
    pfam.taxaInfo.each do |taxon, v|
      out_fh.puts [taxon, v.size].join("\t")
    end
    out_fh.close
    count += 1
    processbar(count, pfams.size)
  end
  puts
end


def create_new_pfams(pfams, type)
  new_pfams = Hash.new
  new2oldPfam = Hash.new
  pfams.each_pair do |pfam_name, v|
    new_pfams[pfam_name.famModify] = v
    new2oldPfam[pfam_name.famModify] = pfam_name
  end
  return([new_pfams, new2oldPfam])
end


def getSpeciesIncluded(infiles)
  species_included = Hash.new
  infiles.each do |infile|
    species_included[getCorename(infile, true)] = ''
  end
  return(species_included)
end


########################################################
if __FILE__ == $0; then
  opts = GetoptLong.new(
    ['--indir', GetoptLong::REQUIRED_ARGUMENT],
    ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
    ['--species_list', GetoptLong::REQUIRED_ARGUMENT],
    ['--seq_indir', GetoptLong::REQUIRED_ARGUMENT],
    ['--seq_suffix', GetoptLong::REQUIRED_ARGUMENT],
    ['--seq_outdir', GetoptLong::REQUIRED_ARGUMENT],
    ['--output_seq', GetoptLong::NO_ARGUMENT],
    ['--fam_include_list', '--fam_list', GetoptLong::REQUIRED_ARGUMENT],
    ['--type', GetoptLong::REQUIRED_ARGUMENT],
    ['-e', '--evalue', GetoptLong::REQUIRED_ARGUMENT],
    ['--cdd_type', '--db_type', GetoptLong::REQUIRED_ARGUMENT],
    ['--bitscore', GetoptLong::NO_ARGUMENT],
    ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
    ['--locus_outdir', GetoptLong::REQUIRED_ARGUMENT],
    ['--super_fam', GetoptLong::NO_ARGUMENT],
    ['--force', GetoptLong::NO_ARGUMENT],
    ['--tolerate', GetoptLong::NO_ARGUMENT],
    ['--all', GetoptLong::NO_ARGUMENT],
    ['--gene2ko', GetoptLong::REQUIRED_ARGUMENT],
  )


  opts.each do |opt, value|
    case opt
      when /^--indir$/
        indirs << value.split(',').map{|i|File.expand_path(i)}
      when /^--outdir/
        outdir = value
      when /^--species_list$/
        species_list_file = value
      when /^--seq_indir$/
        seq_indirs << value.split(',').map{|i|File.expand_path(i)}
      when /^--seq_suffix$/
        seq_suffix = value
      when /^--seq_outdir$/
        seq_outdir = value
      when /^--output_seq$/
        is_output_seq = true
      when /^--(fam_include_list|fam_list)$/
        fam_include_list_files << value.split(',')
      when /^--type$/
        type = value
      when /^(-e|--evalue)$/
        evalue_cutoff = value.to_f
      when /^--(cdd_type|db_type)$/
        cdd_type = value.to_sym
      when /^--cpu$/
        cpu = value.to_i
      when /^--bitscore$/
        is_bitscore = true
      when /^--locus_outdir$/
        locus_outdir = value
      when /^--super_fam$/
        is_super_fam = true
      when /^--force$/
        is_force = true
      when '--tolerate'
        is_tolerate = true
      when /^--all$/
        is_best = false
      when /^--gene2ko$/
        $KEGG_GENE_2_ORTHOs = value.split(',')
        #$KEGG_GENE_2_ORTHOs.flatten!
    end
  end

  if is_output_seq and fam_include_list_files.empty?
    raise "--fam_list has to be specified if you want to output seqs! Exiting ......"
  end

  fam_include_list_files.flatten!


  ########################################################
  if indirs.empty?
    STDERR.puts "--indir has to be specified! Exiting ......" or exit 1
  end



  ########################################################
  type = getTypeBasedOnCddType(type, cdd_type)

  begin
    mkdir_with_force(outdir, is_force, is_tolerate)
  rescue
    STDERR.puts "--outdir has to be specified! Exiting ......"
  end

  if is_output_seq
    seq_outdir = File.join(outdir, 'seq')
    mkdir_with_force(seq_outdir, is_force)
  end

  if not seq_indirs.empty?
    begin
      mkdir_with_force(seq_outdir, is_force)
    rescue
      STDERR.puts "--seq_outdir has to be specified if seq_indir is specified. Exiting ......"
      exit 1
    end
  end


  ########################################################
  mkdir_with_force(locus_outdir, is_force, is_tolerate) unless locus_outdir.nil?

  taxa_included = read_list(species_list_file) if not species_list_file.nil?
  taxa_included = $MODEL_ORGNS unless locus_outdir.nil?

  fam2super = getSuperFamRela($CDD_SUPER_RELA_FILE) if is_super_fam

  if type == 'kegg'
    puts "reading fam ......"
    $KEGG_GENE_2_ORTHOs.each do |kegg_gene_2_ortho|
      # here kegg_gene_2_ortho is a file
      gene2ortho.merge!(getGene2Ortho(kegg_gene_2_ortho))
    end
  else
    cdd_objs = read_cdd_tbl($CDD_TBL)
    cdd_objs = getCddBitScore(cdd_objs, $CDD_BITSCORE_FILE)
    cdd_objs = filterCddType(cdd_objs, cdd_type) unless cdd_type =~ /pfam|hmmsearch/i
  end

  indirs.flatten!.each do |indir|
    infiles << read_infiles(indir)
  end
  infiles.flatten!

  if not seq_indirs.empty?
    puts "reading seq indir ......"
    seq_indirs.flatten!.each do |seq_indir|
      seq_infiles << read_infiles(seq_indir)
    end
    seq_infiles.flatten!

    #species_included = getSpeciesIncluded(infiles)
    seq_objs = read_seq_from_dir(seq_infiles, seq_suffix, taxa_included, cpu, {})

  end

  if not taxa_included.empty?
    infiles.delete_if{|i|not taxa_included.include?(getCorename(File.basename(i)))}
  end

  if not fam_include_list_files.empty?
    fam_include_list_files.each do |file|
      fams_included << read_list(file).keys
    end
  end
  fams_included.flatten!
  #fams_included << read_list(fam_include_list_file).keys if not fam_include_list_files.empty?

  ########################################################
  pfams = runInParallel(type, infiles, cpu, evalue_cutoff, cdd_objs, gene2ortho, seq_objs, fams_included, fam2super, seq_outdir, locus_outdir, taxa_included, is_bitscore, is_best)

  puts pfams.size

  outputResult(pfams, type, outdir)
end


