#! /usr/bin/env ruby


#################################################
require 'getoptlong'
require 'bio'

require 'Dir'


#################################################
infile = nil
outdir = nil
outdir2 = nil
is_force = false
is_tolerate = false


#################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir2', GetoptLong::REQUIRED_ARGUMENT],
  ['--force', GetoptLong::NO_ARGUMENT],
  ['--tolerate', GetoptLong::NO_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when /^-i$/
      infile = value
    when /^--outdir$/
      outdir = value
    when /^--outdir2$/
      outdir2 = value
    when /^--force$/
      is_force = true
    when /^--tolerate$/
      is_tolerate = true
  end
end


#################################################
mkdir_with_force(outdir, is_force, is_tolerate)
mkdir_with_force(outdir2, is_force, is_tolerate)


#################################################
in_fh = File.open(infile, 'r')
in_fh.each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  gene_name, gene_id, gene_accn = line_arr
  outfile = File.join(outdir, gene_id+'.fas')
  if not File.exists?(outfile)
    `esearch -db protein -query "#{gene_accn}" | efetch -format fasta > #{outfile}`
  else
    next if File.size(outfile) == 0
    outfile2 = File.join(outdir2, gene_id+'.fas')
    out_fh = File.open(outfile2, 'w')
    bio_fh = Bio::FlatFile.open(outfile)
    bio_fh.each_entry do |f|
      title = [gene_name, gene_id].join('|')
      out_fh.puts '>'+title
      out_fh.puts f.seq
    end
    bio_fh.close
    out_fh.close
  end
end
in_fh.close


