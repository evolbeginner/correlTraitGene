#! /usr/bin/env ruby


#############################################################
require 'getoptlong'

require 'Dir'


#############################################################
trait_dir = nil
gene_dir = nil
outdir = nil
null_min = nil
null_max = nil
is_force = false
is_title = false
is_bin = false


trait_files = Array.new
traitInfos = Array.new
geneInfoFiles = Array.new
geneInfos = Array.new


#############################################################
class TraitInfo
  attr_accessor :file, :is_bin, :basename, :traitValues
  def initialize(file, is_bin)
    @file = file
    @is_bin = is_bin
    getBasename()
    readTraitFile()
  end

  private
  def getBasename()
    @basename = File.basename(@file)
  end
  def readTraitFile()
    @traitValues = Hash.new
    in_fh = File.open(@file, 'r')
    in_fh.each_line do |line|
      species, value = line.chomp!.split("\t")
      #@traitValues[species] = (value.to_f == value.to_i) ? value.to_i : value.to_f
      @traitValues[species] = value
      if @is_bin
        if @traitValues[species] == '-' or @traitValues[species] == '0'
          ;
        else
          @traitValues[species] = '1'
        end
      end
    end
    in_fh.close
  end
end


#############################################################
def output_traitInfo(traitInfos, geneInfos, outdir, null_min, null_max, is_title)
  traitInfos.each do |traitInfo|
    outdir2 = File.join(outdir, traitInfo.basename)
    mkdir_with_force(outdir2, false, false)
    geneInfos.each do |geneInfo|
      (traitInfo.traitValues.keys - geneInfo.traitValues.keys).each do |species|
        geneInfo.traitValues[species] = '0'
      end
      if not null_max.nil? or not null_min.nil?
        count = (geneInfo.traitValues.keys & traitInfo.traitValues.keys).map{|i|geneInfo.traitValues[i]}.count{|i|i=='0'}
        next if count < null_min unless null_min.nil?
        next if count > null_max unless null_max.nil?
      end
      outfile = File.join(outdir2, geneInfo.basename)
      out_fh = File.open(outfile, 'w')
      out_fh.puts [nil, 'trait', 'predictor'].join("\t") if is_title
      traitInfo.traitValues.each_pair do |species, value|
        out_fh.puts [species, value, geneInfo.traitValues[species]].join("\t")
      end
      out_fh.close
    end
  end
end


#############################################################
if __FILE__ == $0; then

  opts = GetoptLong.new(
    ['--trait_dir', GetoptLong::REQUIRED_ARGUMENT],
    ['--gene_dir', GetoptLong::REQUIRED_ARGUMENT],
    ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
    ['--null_min', GetoptLong::REQUIRED_ARGUMENT],
    ['--null_max', GetoptLong::REQUIRED_ARGUMENT],
    ['--force', GetoptLong::NO_ARGUMENT],
    ['--title', GetoptLong::NO_ARGUMENT],
    ['--bin', GetoptLong::NO_ARGUMENT],
  )


  opts.each do |opt, value|
    case opt
      when /^--trait_dir$/
        trait_dir = value
      when /^--gene_dir$/
        gene_dir = value
      when /^--outdir$/
        outdir = value
      when /^--null_min$/
        null_min = value.to_f
      when /^--null_max$/
        null_max = value.to_f
      when /^--force$/
        is_force = true
      when /^--title$/
        is_title = true
      when /^--bin$/
        is_bin = true
    end
  end



  #############################################################
  mkdir_with_force(outdir, is_force)

  puts "Reading trait files ......"
  trait_files = read_infiles(trait_dir)
  trait_files.select!{|i| File.file?(i) }
  trait_files.map{|i| traitInfos << TraitInfo.new(i, is_bin) }

  puts "Reading gene files ......"
  geneInfoFiles = read_infiles(gene_dir)
  geneInfoFiles.map{|i| geneInfos << TraitInfo.new(i, is_bin) }


  puts "Outputting ......"
  output_traitInfo(traitInfos, geneInfos, outdir, null_min, null_max, is_title)

end


