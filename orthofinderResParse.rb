#! /usr/bin/env ruby


########################################################################
require 'getoptlong'

require 'Dir'
require 'chang_yong'


########################################################################
infile = nil
outfile = nil
include_list_file = nil
is_force = false
is_tolerate = false

taxa = Array.new
famInfo = Hash.new{|h,k|h[k]={}}


########################################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['-o', GetoptLong::REQUIRED_ARGUMENT],
  ['--include_list', GetoptLong::REQUIRED_ARGUMENT],
  ['--force', GetoptLong::NO_ARGUMENT],
  ['--tolerate', GetoptLong::NO_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
    when '-o'
      outfile = value
    when '--include_list'
      include_list_file = value
    when '--force'
      is_force = true
    when '--tolerate'
      is_tolerate = true
  end
end


########################################################################
taxa_included = read_list_file(include_list_file) if not include_list_file.nil?


########################################################################
in_fh = File.open(infile, 'r')
in_fh.each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  if $. == 1
    line_arr[1,line_arr.size-1].each do |i|
      taxa << i
    end
  else
    fam = line_arr[0]
    line_arr[1, line_arr.size-1].each_with_index do |ele, index|
      taxon = taxa[index]
      next if taxon == 'Total'
      next if ele == ''
      famInfo[fam][taxon] = ele.gsub(' ', '')
    end
  end
end
in_fh.close


########################################################################
out_fh = File.open(outfile, 'w')

famInfo.each_pair do |fam, v|
  items = Array.new
  items << fam
  v.each_pair do |taxon, info|
    next if not taxa_included.include?(taxon) if not include_list_file.nil?
    items << [taxon, info].join(':')
  end
  out_fh.puts items.join("\t")
end
out_fh.close


