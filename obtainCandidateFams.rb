#! /usr/bin/env ruby


#################################################################
DIR = File.dirname($0)

$: << File.join(DIR, 'lib')


#################################################################
require 'getoptlong'
require 'parallel'

require 'gainAndLoss'


#################################################################
$ANNOT_FILE = File.expand_path("~/resource/db/cdd/cddid_all.tbl")
$ANNOT_KEGG_FILE = File.expand_path("~/resource/db/kegg/parse/koid_all.tbl")
$ANNOT_RAST_FILE = File.expand_path("~/project/Rhizobiales/results/RAST/results/orthofinder/final.out")
$PNAS_LIST_FILE = File.expand_path("~/project/Rhizobiales/scripts/correlTraitGene/special_prot/pnas_561/pnas.ori.list")
$PNAS_INDIR = File.expand_path("~/project/Rhizobiales/scripts/correlTraitGene/special_prot/pnas_561/analysis")

$GENUS_ORDER = %w[Bradyrhizobium Mesorhizobium Rhizobium Sinorhizobium Bartonella Brucella]
$MODEL_ORGNS = %w[Bradyrhizobium_japonicum_USDA_110 Mesorhizobium_loti_MAFF303099 Rhizobium_etli_CFN_42 Sinorhizobium_meliloti_1021]


#################################################################
infile = nil
lifestyle_infile = nil
type = nil
lr_range = 0..10000
rate_diff_range = -10000..100000
rate_res_range = -10000..10000
is_output_gene_to_fam = false
p_max = 1
q_max = 1
b_max = 1
p_res_max = 1
q_res_max = 1
b_res_max = 1
b_min = 0
cpu = 1


#################################################################
class String
  def famModify
    str = self
    if self =~ /^(PF[^.]+).+/
      str = sub($&, $1).downcase.sub('pf', 'pfam')
    end
    return(str)
  end
end


class FAM
  attr_accessor :func_abbr, :func
  def initialize
    @func_abbr = ''
    @func = ''
  end
end


class BayesTraits
  attr_accessor :fam, :lr, :pvalue, :q, :b, :no_trait0, :prop_trait0, :no_trait1, :prop_trait1, :rate_diff, :rate_res1, :rate_res2, :p_res1, :p_res2, :q_res1, :q_res2, :b_res1, :b_res2
  def initialize(a)
    @fam, @lr, @pvalue, @q, @b, @no_trait0, @prop_trait0, @no_trait1, @prop_trait1, @rate_diff = a
    #@fam, @lr, @pvalue, @q, @b, @no_trait0, @prop_trait0, @no_trait1, @prop_trait1, @rate_res1, @rate_res2, @p_res1, @p_res2, @q_res1, @q_res2, @b_res1, @b_res2 = a
  end
end


#################################################################
def read_pnas_list(list_file)
  locus2gene = Hash.new
  in_fh = File.open(list_file, 'r')
  in_fh.each_line do |line|
    line.chomp!
    gene, locus = line.split("\t")
    locus2gene[locus] = gene
  end
  in_fh.close
  return(locus2gene)
end


def read_pnas_indir(indir, locus2gene)
  fam2gene = Hash.new{|h,k|h[k]=[]}
  gene2fam = Hash.new{|h1,k1|h1[k1] = Hash.new{|h2,k2|h2[k2]=[]} }
  categories = %w[TIGRFAM COG pfam kegg PRK]
  categories.each do |category|
    indir2 = File.join(indir, category)
    Dir.glob("#{indir2}/*").each do |infile|
      fam = File.basename(infile).famModify
      in_fh = File.open(infile, 'r')
      in_fh.each_line do |line|
        line.chomp!
        locus = line.split("\t")[0]
        fam2gene[fam] << locus2gene[locus]
        gene2fam[locus2gene[locus]][category] << fam
        #STDERR.puts fam if fam =~ /^K0/
      end
      in_fh.close
    end
  end
  return([fam2gene, gene2fam])
end


def readLifestyleInfile(infile)
  candidateTaxa = Hash.new
  regexp2genus = {'Bradyrhizobium'=>'Bradyrhizobium', 'Rhizobium'=>'Rhizobium', 'Sinorhizobium'=>'Sinorhizobium', 'Ensifer'=>'Sinorhizobium', 'Agrobacterium'=>'Rhizobium', 'Mesorhizobium'=>'Mesorhizobium', 'Bartonella'=>'Bartonella', 'Brucella'=>'Brucella'}
  candidateGenera = Hash.new{|h,k|h[k]=[]}
  regexp2genus.values.map{|genus| candidateGenera[genus] = [] }

  r = Regexp.new("^Bradyrhizobium|Rhizobium|Sinorhizobium|Ensifer|Agrobacterium|Mesorhizobium|Brucella|Bartonella")
  in_fh = File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    taxon, trait = line.split("\t")
    next if trait != '1'
    next if r !~ taxon
    genus = regexp2genus[$&]
    candidateTaxa[taxon] = genus
    candidateGenera[genus] << taxon
  end
  in_fh.close
  return([candidateTaxa, candidateGenera])
end


def getGenesPresentInCandidateTaxa(infile, candidateTaxa, candidateGenera, cpu)
  countGenus = Hash.new
  indir = File.expand_path(File.join(File.dirname(infile), '../../../'))
  infiles = Array.new
  Dir.glob("#{indir}/*_res/*") do |i|
    infiles << i
  end

  Parallel.map(infiles, in_threads:cpu) do |i|
    b = File.basename(i)
    countGenus[b] = {}
    in_fh = File.open(i, 'r')
    in_fh.each_line do |line|
      line.chomp!
      taxon, pre_ab = line.split("\t").values_at(0,2)
      next if not candidateTaxa.include?(taxon)
      genus = candidateTaxa[taxon]
      countGenus[b][genus] = 0 if not countGenus[b].include?(genus)
      if pre_ab != '0'
        countGenus[b][genus] += 1
      end
      ($GENUS_ORDER - countGenus[b].keys).map{|genus| countGenus[b][genus] = 0}
    end
    in_fh.close
  end

  propGenus = Hash.new{|h,k|h[k]=Hash.new(0)}
  countGenus.each_pair do |fam0, v|
    fam = fam0.famModify
    candidateGenera.keys.map{|genus|propGenus[fam][genus] = 0}
    #v.each_pair do |genus, count|
    $GENUS_ORDER.each do |genus|
      count = v.include?(genus) ? v[genus].to_f : 0
      propGenus[fam][genus] = (count.to_f/candidateGenera[genus].size).round(2)
    end
  end

  return(propGenus)
end


def readAnnotFile(infile)
  famInfo = Hash.new
  in_fh = File.open(infile, 'r')
  in_fh.each_line do |line|
    #214330	CHL00001	rpoB	RNA polymerase beta subunit	1070
    line.chomp!
    line_arr = line.split("\t")
    fam, func_abbr, func = line_arr[1,3]
    famInfo[fam] = FAM.new()
    famInfo[fam].func_abbr = func_abbr
    famInfo[fam].func = func
  end
  in_fh.close
  return(famInfo)
end


def readInfile(infile, lr_range, rate_res_range, rate_diff_range, p_max, q_max, b_max, p_res_max, q_res_max, b_res_max, b_min)
  candidateInfo = Hash.new
  in_fh = File.open(infile, 'r')
  in_fh.each_line do |line|
    #Family  LR  pvalue  qvalue  bonferroni  no_trait0 prop_trait0 no_trait1 prop_trait1 LR_res1 LR_res2 p_res1  p_res2  q_res1  q_res2  bonferroni_res1 bonferroni_res2
    next if $. == 1
    line.chomp!
    line_arr = line.split("\t")
    fam = line_arr[0]
    lr, pvalue, q, b, no_trait0, prop_trait0, no_trait1, prop_trait1, rate_diff = line_arr[1, 9].map{|i|i.to_f}
    #lr, pvalue, q, b, no_trait0, prop_trait0, no_trait1, prop_trait1, rate_diff, rate_res1, rate_res2, p_res1, p_res2, q_res1, q_res2, b_res1, b_res2 = line_arr[1, line_arr.size-1].map{|i|i.to_f}
    next unless lr_range.include?(lr)
    next if pvalue > p_max
    next if q > q_max
    next if b < b_min
    next if b > b_max
    next unless rate_diff_range.include?(rate_diff)
    #next if not (rate_res_range.include?(rate_res1) or rate_res_range.include?(rate_res2))
    #next if p_res1 > p_res_max and p_res2 > p_res_max
    #next if q_res1 > q_res_max and q_res2 > q_res_max
    #next if b_res1 > b_res_max and b_res2 > b_res_max
    candidateInfo[fam] = BayesTraits.new([fam, lr, pvalue, q, b, no_trait0, prop_trait0, no_trait1, prop_trait1, rate_diff])
    #candidateInfo[fam] = BayesTraits.new([fam, lr, pvalue, q, b, no_trait0, prop_trait0, no_trait1, prop_trait1, rate_res1, rate_res2, p_res1, p_res2, q_res1, q_res2, b_res1, b_res2])
  end
  in_fh.close
  return(candidateInfo)
end


def output_gene2fam(locus2gene, gene2fam)
  categories = %w[TIGRFAM COG pfam kegg PRK]
  categories = %w[kegg COG pfam TIGRFAM PRK]
  locus2gene.each_pair do |locus, gene|
    if gene2fam.include?(gene)
      v = gene2fam[gene]
      categories.each do |category|
        print gene+"\t"
        if v.include?(category)
          fams = v[category]
          puts fams.flatten.join("\t")
          break
        end
      end
    else
      puts gene
    end
  end
end



#################################################################
if __FILE__ == $0
  opts = GetoptLong.new(
    ['-i', GetoptLong::REQUIRED_ARGUMENT],
    ['-l', GetoptLong::REQUIRED_ARGUMENT],
    ['--type', GetoptLong::REQUIRED_ARGUMENT],
    ['--lr', GetoptLong::REQUIRED_ARGUMENT],
    ['--rate_res', GetoptLong::REQUIRED_ARGUMENT],
    ['--rate_diff', GetoptLong::REQUIRED_ARGUMENT],
    ['--output_gene2fam', '--output_gene_to_fam', GetoptLong::NO_ARGUMENT],
    ['--p_max', GetoptLong::REQUIRED_ARGUMENT],
    ['--q_max', GetoptLong::REQUIRED_ARGUMENT],
    ['--b_max', GetoptLong::REQUIRED_ARGUMENT],
    ['--p_res_max', GetoptLong::REQUIRED_ARGUMENT],
    ['--q_res_max', GetoptLong::REQUIRED_ARGUMENT],
    ['--b_res_max', GetoptLong::REQUIRED_ARGUMENT],
    ['--b_min', GetoptLong::REQUIRED_ARGUMENT],
    ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
  )


  opts.each do |opt, value|
    case opt
      when /^-i$/
        infile = value
      when /^-l$/
        lifestyle_infile = value
      when /^--type$/
        type = value
      when /^--lr$/
        a = value.split(/[,]/).map{|i|i.to_i}
        lr_range = Range.new(a[0], a[1])
      when '--rate_diff'
        a = value.split(/[,]/).map{|i|i.to_i}
        rate_diff_range = Range.new(a[0], a[1])
      when /^--rate_res$/
        a = value.split(/[,]/).map{|i|i.to_i}
        rate_res_range = Range.new(a[0], a[1])
      when /^(--output_gene2fam|--output_gene_to_fam)$/
        is_output_gene_to_fam = true
      when /^--p_max$/
        p_max = value.to_f
      when /^--q_max$/
        q_max = value.to_f
      when /^--b_max$/
        b_max = value.to_f
      when /^--p_res_max$/
        p_res_max = value.to_f
      when /^--q_res_max$/
        q_res_max = value.to_f
      when /^--b_res_max$/
        b_res_max = value.to_f
      when '--b_min'
        b_min = value.to_f
      when /^--cpu$/
        cpu = value.to_i
    end
  end


  #################################################################
  locus2gene = read_pnas_list($PNAS_LIST_FILE)

  fam2gene, gene2fam = read_pnas_indir($PNAS_INDIR, locus2gene)

  output_gene2fam(locus2gene, gene2fam) and exit if is_output_gene_to_fam

  candidateTaxa, candidateGenera = readLifestyleInfile(lifestyle_infile)

  propGenus = getGenesPresentInCandidateTaxa(infile, candidateTaxa, candidateGenera, cpu)

  famInfo = readAnnotFile($ANNOT_FILE)
  famInfo.merge!(readAnnotFile($ANNOT_KEGG_FILE))
  famInfo.merge!(readAnnotRastFile($ANNOT_RAST_FILE))

  candidateInfo = readInfile(infile, lr_range, rate_res_range, rate_diff_range, p_max, q_max, b_max, p_res_max, q_res_max, b_res_max, b_min)

  puts ['Family', 'pnas_gene', 'function_abbr', 'function', 'pvalue', 'qvalue', 'bonferroni', 'rate_diff', 'prop_trait0', 'prop_trait1', "Brady", "Meso", "Rhizo", "Sino", 'Bartonella', 'Brucella', 'rate_res1', 'rate_res2', 'p_res1', 'p_res2', 'q_res1', 'q_res2', 'b_res1', 'b_res2'].join("\t")
  candidateInfo.sort_by{|fam0,v|v.b}.to_h.each_pair do |fam0, v|
    fam = fam0.famModify
    if not famInfo.include?(fam)
      famInfo[fam] = FAM.new
      if fam =~ /^OG00/
        famInfo[fam].func = ['', '', '', '', '', '']
      end
    end
    puts [fam, fam2gene[fam].join(' '), famInfo[fam].func_abbr, famInfo[fam].func, v.pvalue, v.q, v.b, v.rate_diff, v.prop_trait0, v.prop_trait1, propGenus[fam].sort_by{|i|$GENUS_ORDER.index(i)}.to_h.values, v.rate_res1, v.rate_res2, v.p_res1, v.p_res2, v.q_res1, v.q_res2, v.b_res1, v.b_res2].flatten.join("\t")
  end

end


