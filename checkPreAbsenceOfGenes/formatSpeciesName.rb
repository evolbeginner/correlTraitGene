#! /usr/bin/env ruby


##############################################
# awk -F"\t" 'BEGIN{OFS="\t"}{if($3~/soil/){if($4=="soil"){a="true"}else{a="false"}print $2,a}}' ~/project/Rhizobiales/selection/BioSample/all_info.list > checkPreAbsenceOfGenes/soil.list


##############################################
DIR = File.dirname($0)
$: << File.join(DIR, 'lib')


##############################################
require 'getoptlong'

require 'process_species_name'


##############################################
infile = nil


##############################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when /^-i$/
      infile = value
  end
end


##############################################
in_fh = File.open(infile, 'r')
in_fh.each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  species = line_arr[0]
  species.process_species_name!
  puts line_arr.join("\t")
end
in_fh.close


