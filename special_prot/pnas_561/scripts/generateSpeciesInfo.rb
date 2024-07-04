#! /usr/bin/env ruby


#################################################
require 'getoptlong'


#################################################
infile = nil


#################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when /^-i$/
      infile = value
  end
end


#################################################
in_fh = File.open(infile, 'r')
in_fh.each_line do |line|
  #phaE2	SMC00055	Sinorhizobium meliloti 1021,true	Sinorhizobium meliloti 2011,false
  line.chomp!
  line_arr = line.split("\t")
  gene_name, gene_id = line_arr[0,2]
  species_info_str_arr = line_arr[2, line_arr.size-2]
  
  print [gene_name, gene_id].join("\t") + "\t"
  if not species_info_str_arr.any?{|i|i.split(',')[1] == 'true'}
    puts species_info_str_arr.map{|i|i.split(',')[0]}.join("\t")
  else
    puts species_info_str_arr.select{|i|i.split(',')[1] == 'true'}.map{|i|i.split(',')[0]}.join("\t")
  end
end


