#! /usr/bin/env ruby


################################################################
require 'getoptlong'

require 'chang_yong'


################################################################
ref_file = nil
target_file = nil
field1 = 1
field2s = [1, 2]


################################################################
def read_list_2_array(infile, field)
  items = Array.new
  in_fh = File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    item = line_arr[field-1]
    items << item
  end
  in_fh.close
  return(items)
end


def read_list_2_hash(infile, fields)
  item_sh = Hash.new
  in_fh = File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    item = line_arr[fields[0]-1] #SSW
    item_sh[item] = fields.map{|i| line_arr[i-1] }.join("\t")
  end
  in_fh.close
  return(item_sh)
end


################################################################
opts = GetoptLong.new(
  ['--ref', GetoptLong::REQUIRED_ARGUMENT],
  ['--target', GetoptLong::REQUIRED_ARGUMENT],
  ['--f1', GetoptLong::REQUIRED_ARGUMENT],
  ['--f2', GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '--ref'
      ref_file = value
    when '--target'
      target_file = value
    when '--f1'
      field1 = value.to_i
    when '--f2'
      field2s = value.split(',').map{|i|i.to_i}
  end
end


################################################################
ref_items = read_list_2_array(ref_file, field1)

target_item_sh = read_list_2_hash(target_file, field2s)

ref_items.each do |item|
  puts [item, target_item_sh[item]].join("\t")
end


