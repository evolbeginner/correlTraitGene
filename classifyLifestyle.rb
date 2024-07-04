#! /usr/bin/env ruby


require 'getoptlong'

require 'Dir'


####################################################
infile = nil
outdir = nil
is_force = false


lifestyleInfo = Hash.new
binLifestyle = Hash.new


####################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--force', GetoptLong::NO_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when /^-i$/
      infile = value
    when /^--outdir$/
      outdir = value
    when /^--force$/
      is_force = value
  end
end


mkdir_with_force(outdir, is_force)


####################################################
in_fh = infile == '-' ? STDIN : File.open(infile, 'r')
in_fh.each_line do |line|
  line.chomp!
  species, lifestyle = line.split("\t").values_at(2,3)
  lifestyleInfo[species] = lifestyle
end
in_fh.close if infile != '-'


####################################################
lifestyleInfo.each_pair do |species, lifestyle|
  case lifestyle
    when 'F'
      binLifestyle[species] = [1,0,0,0,1,0,1]
    when 'V'
      binLifestyle[species] = [0,0,1,0,1,0,0]
    when 'R'
      binLifestyle[species] = [0,0,0,1,1,0,'-']
    when 'VR'
      binLifestyle[species] = [0,0,'-','-',1,0,'-']
    when 'M'
      binLifestyle[species] = [0,1,0,0,0,0,'-']
    when 'free living'
      binLifestyle[species] = [0,0,0,0,0,1,'-']
    else
      binLifestyle[species] = %w[- - - - - - -]
  end
end


####################################################
out_fhs = Array.new
outfiles = %W[F-NF M-NM V-NV R-NR P-NP f-Nf F-V].map{|i|File.join(outdir, i)}
out_fhs = outfiles.map{|i|File.open(i, 'w')}

binLifestyle.each_pair do |species, arr|
  arr.each_with_index do |ele, index|
    out_fhs[index].puts [species, ele].join("\t")
  end
end

out_fhs.map{|i|i.close}


