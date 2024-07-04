#! /usr/bin/env ruby


################################################################################
require 'getoptlong'
require 'parallel'
require 'tempfile'

require 'Dir'


################################################################################
$DIR = File.dirname($0)
$PHYLOGLM = File.join($DIR, 'phyloGLM.R')


################################################################################
indir = nil
outfile = nil
outdir = nil
species_tree_file = nil
isBin = nil
cpu = 1
is_force = false


genes = Array.new


################################################################################
class Phyloglm
  attr_accessor :pvalue, :cov, :fdr
  def initialize(gene)
    @gene = gene
    @pvalue = nil
  end
end


class TraitInfo
  attr_accessor :trait, :species, :geneNums
  def initialize(trait)
    @species = Array.new
    @geneNums = Array.new
  end
end


################################################################################
def check_requirements()
  %w[fdr_correction0.py].each do |prog|
    `which #{prog} 2>/dev/null`
    if $?.exitstatus != 0
      raise "Fatal error!\n#{prog} was not installed! Exiting ......"
    end
  end
end


def getAllGenes(indir)
  genes = Array.new
  Dir.foreach(indir) do |b|
    next if b =~ /^\./
    genes << b
  end
  return(genes)
end


def parsePhyloglmOutfile(phyloglmFile, gene)
  phyloglm = Phyloglm.new(gene)
  in_fh = File.open(phyloglmFile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    if line =~ /^predictor/
      line_arr = line.split(/\s+/)
      if line_arr.size >= 7
        phyloglm.cov = line_arr[1].to_f
        phyloglm.pvalue = line_arr[6].to_f
      end
    end
  end
  in_fh.close
  return(phyloglm)
end


def runPhyloglm(genes, indir, outdir, species_tree_file, isBin, cpu)
  phyloglm_info = Hash.new

  results = Parallel.map(genes, in_processes: cpu) do |gene|
    puts gene
    dataframe_file = File.join(indir, gene)
    outfile = File.join(outdir, gene+'.phyloglm')
    `Rscript #{$PHYLOGLM} #{species_tree_file} #{dataframe_file} #{isBin} > #{outfile} 2>/dev/null`
    puts "Rscript #{$PHYLOGLM} #{species_tree_file} #{dataframe_file} #{isBin} > #{outfile} 2>/dev/null"
    phyloglm = parsePhyloglmOutfile(outfile, gene)
    [gene, phyloglm]
  end

  results.each do |arr|
    phyloglm_info[arr[0]] = arr[1]
  end

  return(phyloglm_info)
end


def readDataframeFiles(genes, indir)
  traitInfos = Hash.new

  genes.each do |gene|
    traitInfos[gene] = Hash.new
    infile = File.join(indir, gene)
    in_fh = File.open(infile, 'r')
    in_fh.each_line do |line|
      line.chomp!
      next if line =~ /^\t/
      line_arr = line.split("\t")
      species = line_arr[0]
      trait = line_arr[1]
      if not traitInfos[gene].include?(trait)
        traitInfo = TraitInfo.new(trait)
      else
        traitInfo = traitInfos[gene][trait]
      end
      traitInfo.species << species
      traitInfo.geneNums << line_arr[2].to_i
      traitInfos[gene][trait] = traitInfo
    end
    in_fh.close
  end

  return(traitInfos)
end


def get_fdrs(phyloglm_info)
  pvalues_str = phyloglm_info.sort.to_h.values.map{|i|i.pvalue}.compact.join(',')
  tmp = Tempfile.new("tmp")
  out_fh = File.open(tmp.path, 'w')
  out_fh.puts pvalues_str
  out_fh.close
  #fdrs = `fdr_correction0.py --num #{pvalues_str}`.chomp.split(',').map{|i|i.to_f}
  fdrs = `fdr_correction0.py -i #{tmp.path}`.chomp.split(',').map{|i|i.to_f}
  tmp.unlink

  phyloglm_info.sort.to_h.each_key do |gene|
    phyloglm = phyloglm_info[gene]
    if phyloglm.pvalue.to_s =~ /\d/
      phyloglm_info[gene].fdr = fdrs.shift
    else
      phyloglm_info[gene].fdr = 'NaN'
    end
  end

  return(phyloglm_info)
end



################################################################################
################################################################################
if __FILE__ == $0

  opts = GetoptLong.new(
    ['--indir', GetoptLong::REQUIRED_ARGUMENT],
    ['-o', GetoptLong::REQUIRED_ARGUMENT],
    ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
    ['-s', '--species_tree', GetoptLong::REQUIRED_ARGUMENT],
    ['--bin', GetoptLong::NO_ARGUMENT],
    ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
    ['--force', GetoptLong::NO_ARGUMENT],
  )


  opts.each do |opt, value|
    case opt
      when /^--indir$/
        indir = value
      when /^-o$/
        outfile = value
      when /^--outdir$/
        outdir = value
      when /^(-s|--species_tree)$/
        species_tree_file = value
      when /^--bin$/
        isBin = true
      when /^--cpu$/
        cpu = value.to_i
      when /^--force$/
        is_force = true
    end
  end


  ################################################################################
  raise "outfile has to be given!\n Exiting ......" if outfile.nil?

  check_requirements()

  mkdir_with_force(outdir, is_force)

  genes = getAllGenes(indir)

  traitInfos = readDataframeFiles(genes, indir)


  ################################################################################
  phyloglm_info = runPhyloglm(genes, indir, outdir, species_tree_file, isBin, cpu)

  phyloglm_info = get_fdrs(phyloglm_info)


  ################################################################################
  out_fh = File.open(outfile, 'w')
  phyloglm_info.sort.each do |gene, phyloglm|
    traits = traitInfos[gene].sort.to_h.keys
    trait_str = traits.map{|trait|[traitInfos[gene][trait].geneNums.reduce(:+), (traitInfos[gene][trait].geneNums.reduce(:+)/traitInfos[gene][trait].species.size.to_f).round(2)].join("\t")}.join("\t")
    out_fh.puts [gene, phyloglm.cov, phyloglm.pvalue, phyloglm.fdr, trait_str].join("\t")
  end
  out_fh.close

end


