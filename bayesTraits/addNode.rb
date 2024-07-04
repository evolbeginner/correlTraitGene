#! /usr/bin/env ruby


#########################################################
$: << File.join(Dir.pwd, 'lib')


#########################################################
require 'getoptlong'

require 'bio'

require 'tree'


#########################################################
tree_file = nil
nexus_outfile = nil
tag_outfile = nil
modes = Array.new
isAddNode = false
sct = 0.001


tags_arr = Array.new


#########################################################
def formatTagNode(nodes, tag_name, node_name)
  #AddTag Tag-2 Sheep Goat
  #AddNode Node-2 Tag-2
  rv_arr = Array.new
  node_str = nodes.map{|i|i.name.gsub(' ', '_')}.join("\t")
  rv_arr << ["AddTag", tag_name, node_str].join("\t")
  rv_arr << ["AddNode", node_name, tag_name].join("\t")
  return(rv_arr)
end


def output_tag(tag_outfile, tags_arr, modes, sct)
  out_fh = File.open(tag_outfile, 'w')
  modes.each do |mode|
    out_fh.puts mode
  end
  out_fh.puts

  out_fh.puts "SCT #{sct}"

  tags_arr.each do |arr|
    out_fh.puts arr.join("\n")
  end

  out_fh.puts
  out_fh.puts 'run'
  out_fh.close
end


#########################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--nexus', GetoptLong::REQUIRED_ARGUMENT],
  ['--addNode', GetoptLong::NO_ARGUMENT],
  ['--tag', GetoptLong::REQUIRED_ARGUMENT],
  ['--mode', GetoptLong::REQUIRED_ARGUMENT],
  ['--sct', '--SCT', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when /^-i$/
      tree_file = value
    when /^--nexus$/
      nexus_outfile = value
    when /^--addNode$/
      isAddNode = true
    when /^--tag$/
      tag_outfile = value
    when /^--mode$/
      modes = value.split(',')
    when /^--sct$/i
      sct = value.to_f
  end
end


#########################################################
if nexus_outfile.nil?
  puts "Fatal error! Nexus_outfile is not provided!"
  exit 1
end


if not tag_outfile.nil? and not modes.empty?
  isAddNode = true
elsif not tag_outfile.nil? and modes.empty?
  puts "Wrong!"
  exit 1
elsif tag_outfile.nil? and not modes.empty?
  puts "Wrong!"
  exit 1
end


#########################################################
trees = getTreeObjs(tree_file)
tree = trees[0]
tree.options[:bootstrap_style] = :disabled
$stdout = File.open(nexus_outfile, 'w')
tree.outputNexus(true, true)
$stdout = STDOUT


#########################################################
exit if not isAddNode

tree.internal_nodes.each_with_index do |node, index|
  nodes = tree.tips(node)
  tag_name = ['Tag', (index+1).to_s].join('-')
  node_name = ['Node', (index+1).to_s].join('-')
  tags_arr << formatTagNode(nodes, tag_name, node_name)
end

output_tag(tag_outfile, tags_arr, modes, sct)


