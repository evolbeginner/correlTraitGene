#! /usr/bin/env ruby


#####################################################
require 'getoptlong'
require 'bio'

require 'util'


#####################################################
indir = nil


#####################################################
opts = GetoptLong.new(
  ['--indir', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '--indir'
      indir = value
  end
end


#####################################################
Dir.foreach(indir) do |b|
  next if b =~ /^\./
  infile = File.join(indir, b)
  c = getCorename(infile)
  in_fh = Bio::FlatFile.open(infile)
  in_fh.each_entry do |f|
    title = f.definition.split(' ')[0]
    puts [title, c].join("\t")
  end
  in_fh.close
end


