#! /usr/bin/env ruby


#####################################################
require 'getoptlong'
require 'parallel'

require 'chang_yong'
require 'Dir'
require 'util'


#####################################################
infiles = Array.new
gene_indir = nil
mins = Array.new
maxs = Array.new
cpu = 1

a = Array.new
gene2num = Hash.new{|h1,k1|h1[k1]=Hash.new{|h2,k2|h2[k2]=0}}
ranges = Array.new


#####################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--i1', GetoptLong::REQUIRED_ARGUMENT],
  ['--i2', GetoptLong::REQUIRED_ARGUMENT],
  ['--minmax', GetoptLong::REQUIRED_ARGUMENT],
  ['--indir', '--gene_indir', GetoptLong::REQUIRED_ARGUMENT],
  ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '-i'
      infiles << infiles.split(',')
    when '--i1'
      infiles << value
    when '--i2'
      infiles << value
    when '--minmax'
      value.split(',').each do |v|
        min, max = v.split('-').map(&:to_f)
        mins << min
        maxs << max
      end
    when '--indir', '--gene_indir'
      gene_indir = value
    when '--cpu'
      cpu = value.to_i
  end
end


infiles.flatten!


#####################################################
mins.zip(maxs).each do |min, max|
  ranges << Range.new(min, max)
end



#####################################################
infiles.each do |infile|
  a << read_list(infile).keys
end


#####################################################
gene_infiles = read_infiles(gene_indir)
gene_infiles.select!{|i|File.file?(i)}


results = Parallel.map(gene_infiles, in_processes: cpu) do |gene_infile|
  c = getCorename(gene_infile)
  gene2num = Hash.new
  gene2num[c] = {0=>0, 1=>0}
  in_fh = File.open(gene_infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    taxon, no = line.split("\t")
    no = no.to_i
    next if no == 0
    if a[0].include?(taxon)
      gene2num[c][0] += 1
    end
    if a[1].include?(taxon)
      gene2num[c][1] += 1
    end
  end
  in_fh.close
  gene2num
end


##################################################
results.each do |h|
  gene2num.merge!(h)
end


##################################################
gene2num.each_pair do |gene, v|
  prop1 = v[0]/a[0].size.to_f
  prop2 = v[1]/a[1].size.to_f
  if ranges[0].include?(prop1) and ranges[1].include?(prop2)
    puts [gene, prop1, prop2].join("\t")
  end
end


