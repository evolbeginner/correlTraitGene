#! /usr/bin/env ruby


################################################################
require 'getoptlong'

require 'Dir'
require 'chang_yong'


################################################################
indir = nil
fam_list_file = nil

fams_included = Hash.new
taxon2fam = Hash.new{|h,k|h[k]={}}


################################################################
def getGenePreAbsenceInfo(infile, taxon2fam)
  b = File.basename(infile)
  in_fh = File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    taxon, v = line.split("\t").values_at(0,2)
    taxon2fam[taxon][b] = v.to_i
  end
  in_fh.close
  return(taxon2fam)
end


################################################################
opts = GetoptLong.new(
  ['--indir', GetoptLong::REQUIRED_ARGUMENT],
  ['--fam_list', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '--indir'
      indir = value
    when '--fam_list'
      fam_list_file = value
  end
end


################################################################
infiles = read_infiles(indir)

fams_included = read_list(fam_list_file)

infiles.each do |infile|
  b = File.basename(infile)
  next unless fams_included.include?(b)
  taxon2fam = getGenePreAbsenceInfo(infile, taxon2fam)
end


################################################################
puts ['species', taxon2fam[taxon2fam.keys[0]].keys].flatten.join("\t")
taxon2fam.each_pair do |taxon, v|
  puts [taxon, v.values].flatten.join("\t")
end


