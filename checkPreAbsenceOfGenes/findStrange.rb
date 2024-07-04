#! /usr/bin/env ruby


############################################################
require 'getoptlong'


############################################################
infile = nil
lifestyles = Array.new
lifestyle_pre_ab = nil
gene_pre_ab = nil


count_info = Hash.new


############################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--lifestyle', GetoptLong::REQUIRED_ARGUMENT],
  ['--lifestyle_pre_ab', GetoptLong::REQUIRED_ARGUMENT],
  ['--gene_pre_ab', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when /^-i$/
      infile = value
    when /^--lifestyle$/
      lifestyles << value.split(',')
    when /^--lifestyle_pre_ab/
      lifestyle_pre_ab = value
    when /^--gene_pre_ab$/
      gene_pre_ab = value
  end
end


lifestyles.flatten!


############################################################
in_fh = File.open(infile, 'r')
in_fh.each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  taxon, lifestyle, is_lifestyle = line_arr.values_at(0,1,2)
  if lifestyles.include?(lifestyle) or is_lifestyle == lifestyle_pre_ab
    count_info[taxon] = line_arr[3, line_arr.size-3].count(gene_pre_ab)
  end
end
in_fh.close


############################################################
count_info.each_pair do |taxon, count|
  puts [taxon, count].join("\t")
end


