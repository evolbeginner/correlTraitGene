#! /usr/bin/env ruby


################################################
require 'getoptlong'
require 'parallel'


################################################
infile = nil
cpu = 1


################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
    when '--cpu'
      cpu = value.to_i
  end
end


################################################
in_fh = File.open(infile, 'r')
lines = in_fh.readlines.map{|i|i.chomp}
in_fh.close


Parallel.map(lines, in_threads: cpu) do |line|
  line_arr = line.split("\t")
  gene, id = line_arr.values_at(0, 1)
  n = `esearch -db gene -query #{id} | efetch -format tabular | grep -i #{id} | grep -i #{gene} | wc -l 2>/dev/null`.chomp.to_i
  if n != 1
    puts [gene, id, n].join("\t")
  end
end


