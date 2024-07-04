#! /usr/bin/env ruby


################################################################
require 'getoptlong'

require 'Dir'
require 'chang_yong'


################################################################
list_files = Array.new
indirs = Array.new
lifestyle_file = nil

infiles = Array.new
target_infiles = Array.new
taxon_info = Hash.new
domains_included = Array.new


################################################################
class Taxon
  attr_accessor :name, :lifestyle, :is_lifestyle, :presences
  def initialize()
    @presences = Hash.new{|h,k|h[k]=[]}
  end
end


################################################################
def read_lifestyle_file(infile)
  lifestyle_info = Hash.new
  in_fh = File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    species, lifestyle = line_arr.values_at(2,3)
    lifestyle_info[species] = lifestyle
  end
  in_fh.close
  return(lifestyle_info)
end


################################################################
opts = GetoptLong.new(
  ['--indir', GetoptLong::REQUIRED_ARGUMENT],
  ['--list', GetoptLong::REQUIRED_ARGUMENT],
  ['--lifestyle_file', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when /^--list$/
      list_files << value.split(',')
    when /^--indir$/
      indirs << value.split(',')
    when /^--lifestyle_file/
      lifestyle_file = value
  end
end


list_files.flatten!
indirs.flatten!


################################################################
lifestyle_info = read_lifestyle_file(lifestyle_file)

list_files.each do |list_file|
  domains_included << read_list(list_file).keys
end


indirs.each do |indir|
  infiles << read_infiles(indir)
end
infiles.flatten!


################################################################
infiles.each do |infile|
  b = File.basename(infile)
  if domains_included.any?{|i|i.include?(b)}
    target_infiles << infile
  end
end


target_infiles.each do |target_infile|
  b = File.basename(target_infile)
  in_fh = File.open(target_infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    taxon_name, is_lifestyle, is_present = line_arr.values_at(0,1,2)
    if taxon_info.include?(taxon_name)
      taxon = taxon_info[taxon_name]
    else
      taxon = Taxon.new
    end
    taxon.name = taxon_name
    taxon.is_lifestyle = is_lifestyle
    taxon.lifestyle = lifestyle_info[taxon_name]
    domains_included.each_with_index do |domains, index|
      taxon.presences[index] << is_present if domains.include?(b)
    end
    taxon_info[taxon.name] = taxon
  end
  in_fh.close
end


################################################################
taxon_info.sort.to_h.each_pair do |taxon_name, taxon|
  next if taxon.lifestyle.nil?
  puts [taxon.name, taxon.lifestyle, taxon.is_lifestyle, taxon.presences.map{|type, arr|arr}].flatten.join("\t")
end


