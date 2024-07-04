#! /usr/bin/env ruby


########################################################################
require 'getoptlong'

require 'Dir'


########################################################################
infile = nil
outdir = nil
is_force = false
is_tolerate = false

taxa = Array.new
famInfo = Hash.new{|h,k|h[k]={}}


########################################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--force', GetoptLong::NO_ARGUMENT],
  ['--tolerate', GetoptLong::NO_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
    when '--outdir'
      outdir = value
    when '--force'
      is_force = true
    when '--tolerate'
      is_tolerate = true
  end
end


########################################################################
mkdir_with_force(outdir, is_force, is_tolerate)


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
    line_arr[1, line_arr.size-1].each_with_index do |num, index|
      taxon = taxa[index]
      next if taxon == 'Total'
      famInfo[fam][taxon] = num.to_i
    end
  end
end
in_fh.close


########################################################################
famInfo.each_pair do |fam, v|
  outfile = File.join(outdir, fam)
  out_fh = File.open(outfile, 'w')

  v.each_pair do |taxon, num|
    out_fh.puts [taxon, num].join("\t")
  end
  out_fh.close
end


