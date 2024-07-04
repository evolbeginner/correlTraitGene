#! /usr/bin/env ruby


############################################
require 'getoptlong'

require 'Dir'
require 'processbar'
require 'chang_yong'


############################################
infile = nil
outfile = nil
outdir = nil
genome_outdir = nil
seq_outdir = nil
species_list_file = nil
is_force = false
is_tolerate = false


############################################
class Gene
  attr_accessor :gene, :gene_id, :taxon, :is_discon, :accn, :coor
end


############################################
def getGeneIds(infile)
  geneIds = Hash.new
  in_fh = File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    gene, gene_id = line.split("\t")
    geneIds[gene_id] = gene
  end
  in_fh.close
  return(geneIds)
end


def inBatchSearch(geneIds, outdir, species_list_file)
  gene_id_2_taxon = Hash.new{|h,k|h[k]=[]}
  count = 0
  total = geneIds.keys.size

  geneIds.each_pair do |gene_id, gene|
    count += 1
    outfile = File.join(outdir, gene_id)
    if not File.exists?(outfile)
      system("esearch -db gene -query #{gene_id} | efetch > #{outfile}")
    end
    output_arr = `cat #{outfile}`.chomp.split(/\n\n/)
    output_arr.each do |item_str|
      items = item_str.split("\n").select{|i|i!=''}
      items[1] =~ /\[(.+)\]/
      taxon = $1
      next if not items[0].split('. ')[1].include?(gene)
      if not items.select{|i|i == 'This record was discontinued.'}.empty?
        is_discon = false
      else
        is_discon = true
      end

      #Annotation:  NC_009937.1 (11440..11892)
      accn = nil
      coor = nil
      items.each do |line|
        if line =~ /^Annotation:.+(\w\w_[^ ]+) \((.+)\)/
          accn = $1
          coor = $2
        end
      end      

      gene_item = Gene.new
      gene_item.gene = gene
      gene_item.gene_id = gene_id
      gene_item.taxon = taxon
      gene_item.is_discon = is_discon
      gene_item.accn = accn
      gene_item.coor = coor
      gene_id_2_taxon[gene_id] << gene_item
    end
    processbar(count, total)
  end
  puts
  return(gene_id_2_taxon)
end


def output_result(gene_id_2_taxon, geneIds, outfile)
  out_fh = File.open(outfile, 'w')
  gene_id_2_taxon.each_pair do |gene_id, taxaInfo|
    gene = geneIds[gene_id]
    out_fh.puts [gene, gene_id, taxaInfo.map{|i|[i.taxon, i.is_discon].join(',')}].flatten.join("\t")
  end
  out_fh.close
end


def output_accn(gene_id_2_taxon, geneIds, species, genome_outdir, seq_outdir)
  gene_id_2_taxon.each_pair do |gene_id, taxaInfo|
    gene = geneIds[gene_id]
    puts [gene, gene_id, taxaInfo.select{|i|species.include?(i.taxon)}.map{|i|[i.accn, i.coor].join('|')}].flatten.join("\t")
    good_gene_items = taxaInfo.select{|i|species.include?(i.taxon)}
    if good_gene_items[0].nil?
      p ["Problem!\t", taxaInfo]
      next
    end
    [good_gene_items[0]].each do |i|
      accn = i.accn
      coor = i.coor
      genome_outfile = File.join(genome_outdir, accn+'.fas')
      if not File.exists?(genome_outfile)
        `esearch -db nucleotide -query \'#{accn}[accn]\' | efetch -format fasta > #{genome_outfile}`
      end
      bed_outfile = 'test.bed'
      generate_bed_file(i, bed_outfile)
      seq_outfile = File.join(seq_outdir, gene_id+'.fas')
      `bedtools getfasta -fi #{genome_outfile} -bed test.bed -fo #{seq_outfile} -name -s`
    end
  end
end


def generate_bed_file(i, bed_outfile)
  out_fh = File.open(bed_outfile, 'w')
  strand = '+'
  if i.coor =~ /complement/i
    strand = '-'
  end
  i.coor =~ /^(\d+)\.\.(\d+)/
  start, stop = $1.to_i, $2.to_i

  if strand == '+'
    out_fh.puts [i.accn, start-1, stop-1, i.gene+'|'+i.gene_id, '1', strand].join("\t")
  else
    out_fh.puts [i.accn, start, stop, i.gene+'|'+i.gene_id, '1', strand].join("\t")
  end
  
  out_fh.close
end


############################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['-o', GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--genome_outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--seq_outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--species_list', GetoptLong::REQUIRED_ARGUMENT],
  ['--force', GetoptLong::NO_ARGUMENT],
  ['--tolerate', GetoptLong::NO_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when /^-i$/
      infile = value
    when /^-o$/
      outfile = value
    when /^--outdir$/
      outdir = value
    when /^--genome_outdir$/
      genome_outdir = value
    when /^--seq_outdir$/
      seq_outdir = value
    when /^--species_list$/
      species_list_file = value
    when /^--force$/
      is_force = true
    when /^--tolerate$/
      is_tolerate = true
  end
end


############################################
species = read_list(species_list_file) if not species_list_file.nil?

mkdir_with_force(outdir, is_force, is_tolerate)

geneIds = getGeneIds(infile)

gene_id_2_taxon = inBatchSearch(geneIds, outdir, species)


if not species_list_file.nil?
  mkdir_with_force(genome_outdir, is_force, is_tolerate)
  mkdir_with_force(seq_outdir, is_force, is_tolerate)
  output_accn(gene_id_2_taxon, geneIds, species, genome_outdir, seq_outdir)
else
  output_result(gene_id_2_taxon, geneIds, outfile)
end


