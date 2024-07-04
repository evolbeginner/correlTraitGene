#! /usr/bin/env ruby


##################################################################
require 'getoptlong'
require 'parallel'


require 'Dir'
require 'chang_yong'


##################################################################
indir = nil
cpu = 1
include_list_file = nil


fam2taxon_info = Hash.new
taxon2fam_info = Hash.new
species_included = Hash.new


##################################################################
def getGeneNum(infile, species_included)
  fam2taxon_info = Hash.new
  taxon2fam_info = Hash.new
  in_fh = File.open(infile, 'r')
  fam_name = File.basename(infile)
  in_fh.each_line do |line|
    line.chomp!
    taxon, fam_num = line.split("\t")
    next unless species_included.include?(taxon) unless species_included.empty?
    fam_num = fam_num.to_i
    fam2taxon_info[fam_name] = Hash.new if not fam2taxon_info.include?(fam_name)
    taxon2fam_info[taxon] = Hash.new if not taxon2fam_info.include?(taxon)
    fam2taxon_info[fam_name][taxon] = fam_num
    taxon2fam_info[taxon][fam_name] = fam_num
  end
  in_fh.close
  return([fam2taxon_info, taxon2fam_info])
end


##################################################################
opts = GetoptLong.new(
  ['--indir', GetoptLong::REQUIRED_ARGUMENT],
  ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
  ['--include_list', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when /^--indir$/
      indir = value
    when /^--cpu$/
      cpu = value.to_i
    when /^--include_list$/
      include_list_file = value
  end
end


##################################################################
infiles = read_infiles(indir)

species_included = read_list(include_list_file)

results = Parallel.map(infiles, in_processes: cpu) do |infile|
  fam2taxon_info, taxon2fam_info = getGeneNum(infile, species_included) if File.file?(infile)
  [fam2taxon_info, taxon2fam_info]
end

results.each do |fam2taxon_info_1, taxon2fam_info_1|
  fam2taxon_info.merge!(fam2taxon_info_1)
  taxon2fam_info.merge!(taxon2fam_info_1)
end


##################################################################
taxa = taxon2fam_info.keys.sort

puts ['', taxa].flatten.join("\t")
fam2taxon_info.each_pair do |fam_name, v|
  copy_num_str = taxa.map{|taxon|a = v.include?(taxon) ? v[taxon] : 0}.join("\t")
  puts [fam_name, copy_num_str].join("\t")
end


