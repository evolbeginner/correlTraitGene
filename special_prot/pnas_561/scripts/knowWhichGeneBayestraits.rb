#! /usr/bin/env ruby


##########################################################
require 'getoptlong'

require 'Dir'


##########################################################
indir = nil
summary_file = nil
list_file = nil

basename_info = Hash.new{|h,k|h[k]=[]}


##########################################################
def get_target_gene(infile)
  func_info = Hash.new
  in_fh = File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    func_id = line_arr[0]
    func_info[func_id] = ''
  end
  in_fh.close
  return(func_info)
end


##########################################################
opts = GetoptLong.new(
  ['--indir', GetoptLong::REQUIRED_ARGUMENT],
  ['--summary', GetoptLong::REQUIRED_ARGUMENT],
  ['--list', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when /^--indir$/
      indir = value
    when /^--summary$/
      summary_file = value
    when /^--list$/
      list_file = value
  end
end


##########################################################
#for i in TIGRFAM/TIGRFAM_res/*; do b=`basename $i`; if grep $b ~/project/Rhizobiales/scripts/correlTraitGene/TIGRFAM/bayestraits/summary/F-NF.txt >/dev/null; then echo -ne "$b\t"; awk '{print $1}' $i; fi done

infiles = read_infiles(indir)
infiles.each do |infile|
  b = File.basename(infile)
  in_fh = File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    gene_id = line.split("\t")[0]
    basename_info[b] << gene_id
  end
  in_fh.close
end

func_info = get_target_gene(summary_file)

(basename_info.keys & func_info.keys).each do |func_id|
  puts [func_id, basename_info[func_id]].flatten.join("\t")
end


