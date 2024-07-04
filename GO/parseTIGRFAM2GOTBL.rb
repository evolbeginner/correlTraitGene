#! /usr/bin/env ruby


#############################################################
require 'getoptlong'


#############################################################
infile = nil


#############################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when /^-i$/
      infile = value
  end
end


#############################################################
in_fh = File.open(infile, 'r')
in_fh.each_line do |line|
  next if line =~ /^!/
  #JCVI_TIGRFAMS:TIGR00001 ribosomal protein L35 > GO:translation ; GO:0006412
  line.chomp!
  tigrName = line.split(' ')[0].split(':')[1]
  if line =~ /GO:(.+) ;/
    goName = $1
  end
  if line =~ /GO:(\S+)$/
    goId= $&
  end
  #puts [tigrName, goId].join("\t")
  #role_id:	100	mainrole:	Central intermediary metabolism
  puts ['role_id:', goId, '', goName].join("\t")
end
in_fh.close


