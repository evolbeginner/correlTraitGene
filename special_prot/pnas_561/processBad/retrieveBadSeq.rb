#! /usr/bin/env ruby


#########################################################
require 'getoptlong'
require 'parallel'

require 'Dir'


#########################################################
infile = nil
cpu = 1
outdir = nil
is_force = true


#########################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--force', GetoptLong::NO_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when /^-i$/
      infile = value
    when /^--cpu$/
      cpu = value.to_i
    when /^--outdir$/
      outdir = value
    when /^--force$/
      is_force = true
  end
end


#########################################################
in_fh = File.open(infile, 'r')
lines = in_fh.readlines.map{|i|i.chomp}
in_fh.close


mkdir_with_force(outdir, is_force)


#########################################################
Parallel.map(lines, in_processes:cpu) do |line|
  line_arr = line.split("\t")
  gene, id, accn = line_arr.values_at(0,1,3)
  #p "esearch -db protein -query #{accn} | efetch -format fasta"
  outfile = File.join(outdir, id+'.protein')
  out_fh = File.open(outfile, 'w')
  seq = `esearch -db protein -query "#{accn}" 2>/dev/null | efetch -format fasta 2>/dev/null`
  p "esearch -db protein -query #{accn} 2>/dev/null | efetch -format fasta 2>/dev/null"
  out_fh.puts seq
  out_fh.close
end


