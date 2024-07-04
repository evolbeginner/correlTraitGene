#! /usr/bin/env ruby


#####################################################
require 'getoptlong'
require 'bio'

require 'Dir'
require 'SSW_bio'
require 'util'


#####################################################
list_infiles = Array.new
indirs = Array.new

seq_objs = Hash.new{|h,k|h[k]={}}
locus2gene = Hash.new
infiles = Array.new


#####################################################
def read_rela_file(infile)
  a2b = Hash.new
  b2a = Hash.new
  in_fh = File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    a, b = line_arr.values_at(0,1)
    a2b[a] = b
    b2a[b] = a
  end
  in_fh.close
  return([a2b, b2a])
end


#####################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--indir', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '-i'
      list_infiles << value.split(',')
    when '--indir'
      indirs << value.split(',')
  end
end


list_infiles.flatten!
indirs.flatten!


#####################################################
list_infiles.each do |list_infile|
  locus2gene.merge! read_rela_file(list_infile)[1]
end

indirs.each do |indir|
  infiles.concat read_infiles(indir)
end

infiles.each do |infile|
  c = getCorename(infile)
  read_seq_file(infile).each_pair do |title, f|
    seq_objs[c] = f
  end
end


#####################################################
seq_objs.each_pair do |c, f|
  accn, annot, taxon = nil, nil, nil
  if f.definition =~ /([^ ]+) ([^\[]+) \[(.+)\]$/
    #WP_012169134.1 hydrogenase [Azorhizobium caulinodans]
    accn, annot, taxon = $1, $2, $3
  elsif f.definition =~ /([^ ]+) (.+) OS=([^=\(\)]+) /
    #>tr|B0ZTC1|B0ZTC1_RHILI Acyl-homoserine-lactone synthase OS=Rhizobium loti OX=381 GN=mrlI2 PE=3 SV=1
    accn, annot, taxon = $1, $2, $3
  elsif f.definition =~ /([^ ]+) (.+)/
    accn, annot = $1, $2
  end
  puts [locus2gene[c], c, accn, annot, taxon].join("\t")
end


