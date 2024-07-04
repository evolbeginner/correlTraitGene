#! /usr/bin/env ruby


##################################################################
$: << File.expand_path("~/project/Rhizobiales/scripts/angst/")


##################################################################
require 'getoptlong'
require 'parallel'

require 'chang_yong'
require 'Array'
require 'enrichedRast-repeated'


##################################################################
$FISHER = 'fisher_chi_test.py'
$FDR_CORR = 'fdr_correction0.py'

TYPES = [:categories, :subcategories, :subsystems]


##################################################################
infile = nil
background = nil
rast_file = nil
cpu = 1


pvalue_info = Hash.new{|h,k|h[k]={}}
fdr_info = Hash.new{|h,k|h[k]={}}


##################################################################
def getBackgroundFams(background)
  items = Hash.new
  if File.ftype(background) == 'file'
    items = read_list(background)
  elsif File.ftype(background) == 'directory'
    Dir.foreach(background).each do |b|
      next if b =~ /^\./
      items[b] = ''
    end
  end
  return(items)
end


def getFams(infile)
  fams = Hash.new
  in_fh = File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    fam = line_arr[0]
    fams[fam] = ''
  end
  in_fh.close
  return(fams)
end


def get_final(fams, rastInfo)
  final = Hash.new{|h,k|h[k]=[]}
  TYPES.each do |type|
    fams.each_key do |fam|
      next if not rastInfo.include?(fam)
      rast = rastInfo[fam]
      final[type] << rast.send(type)
    end
    final[type].flatten!
  end
  return(final)
end


##################################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['-b', GetoptLong::REQUIRED_ARGUMENT],
  ['--rast', GetoptLong::REQUIRED_ARGUMENT],
  ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
    when '-b'
      background = value
    when '--rast'
      rast_file = value
    when '--cpu'
      cpu = value.to_i
  end
end


##################################################################
rastInfo = readRastFile(rast_file)

background_fams = getBackgroundFams(background)

fams = getFams(infile)


##################################################################
final1 = get_final(fams, rastInfo)
final2 = get_final(background_fams, rastInfo)


TYPES.each do |type|
  freqs1 = final1[type].getFreqs
  freqs2 = final2[type].getFreqs
  sum1 = freqs1.values.sum
  sum2 = freqs2.values.sum
  Parallel.map(freqs1, in_threads: cpu) do |annot, n1|
    next if annot == ''
    n2 = freqs2[annot]
    num_str = [n1, n2, sum1, sum2].join(',')
    pvalue = `#{$FISHER} --type fisher --num #{num_str}`.chomp.to_f
    pvalue_info[type][annot] = pvalue
  end

  pvalues = pvalue_info[type].values.select{|i|i!=1}
  fdrs = `#{$FDR_CORR} --num #{pvalues.join(',')}`.chomp.split(',').map{|i|i.to_f}
  pvalue_info[type].select{|k,v|v!=1}.keys.zip(fdrs) do |annot, fdr|
    fdr_info[type][annot] = fdr
  end
end


fdr_info.each_pair do |type, v1|
  freqs1 = final1[type].getFreqs
  freqs2 = final2[type].getFreqs
  sum1 = freqs1.values.sum
  sum2 = freqs2.values.sum
  v1.each_pair do |annot, fdr|
    pv = pvalue_info[type][annot]
    puts [type, annot, freqs1[annot], freqs2[annot], sum1, sum2, pv, fdr].join("\t")
  end
end


