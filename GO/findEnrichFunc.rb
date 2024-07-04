#! /usr/bin/env ruby


##########################################################
require 'getoptlong'
require 'parallel'


##########################################################
#ROLE_LINK_FILE = File.expand_path("~/resource/db/TIGRFAM/TIGRFAMS_ROLE_LINK")
#ROLE_NAME_FILE = File.expand_path("~/resource/db/TIGRFAM/TIGR_ROLE_NAMES")
ROLE_LINK_FILE = File.expand_path("~/resource/db/TIGRFAM/TIGRFAM_GO_ID")
ROLE_NAME_FILE = File.expand_path("~/resource/db/TIGRFAM/TIGRFAM_GO_NAME")

infile = nil
cpu = 1
is_output_role = false

pvalue_info = Hash.new
fdr_info = Hash.new


##########################################################
def get_TIGRFAM_roles(role_link_file, role_name_file)
  tigrName2RoleId = Hash.new
  roleId2TigrNames = Hash.new{|h,k|h[k]=[]}
  in_fh = File.open(role_link_file, 'r')
  in_fh.each_line do |line|
    line.chomp!
    tigrName, roleId = line.split("\t")
    roleId = roleId
    tigrName2RoleId[tigrName] = roleId
    roleId2TigrNames[roleId] << tigrName
  end
  in_fh.close

  roleId2RoleName = Hash.new
  roleName2RoleId = Hash.new
  roleLevelInfo = Hash.new
  in_fh = File.open(role_name_file, 'r')
  in_fh.each_line do |line|
    #role_id:	100	mainrole:	Central intermediary metabolism
    line.chomp!
    roleId, roleLevel, roleName = line.split("\t").values_at(1,2,3)
    roleId = roleId
    roleId2RoleName[roleId] = roleName
    roleName2RoleId[roleName] = roleId
    roleLevelInfo[roleName] = roleLevel
  end
  in_fh.close

  roleName2TigrNames = Hash.new{|h,k|h[k]=[]}
  tigrName2RoleNames = Hash.new{|h,k|h[k]=[]}
  roleId2TigrNames.each_pair do |roleId, tigrNames|
    roleName = roleId2RoleName[roleId]
    next if roleName.nil?
    roleName2TigrNames[roleName] = tigrNames
    tigrNames.each do |tigrName|
      tigrName2RoleNames[tigrName] << roleName
    end
  end

  return([roleName2TigrNames, tigrName2RoleNames, roleLevelInfo])
end


def readInfile(infile, tigrName2RoleNames)
  tigrNamesInput = Hash.new{|h,k|h[k]=[]}
  in_fh = File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    next if $. == 1
    tigrName = line.split("\t")[0]
    roleNames = tigrName2RoleNames[tigrName]
    next if not tigrName2RoleNames.include?(tigrName)
    roleNames.each do |roleName|
      tigrNamesInput[roleName] << tigrName
    end
  end
  in_fh.close
  return(tigrNamesInput)
end


##########################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
  ['--output_role', GetoptLong::NO_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when /^-i$/
      infile = value
    when /^--cpu$/
      cpu = value.to_i
    when /^--output_role$/
      is_output_role = true
  end
end


##########################################################
roleName2TigrNames, tigrName2RoleNames, roleLevelInfo = get_TIGRFAM_roles(ROLE_LINK_FILE, ROLE_NAME_FILE)

tigrNamesInput = readInfile(infile, tigrName2RoleNames)

#total_size = roleName2TigrNames.values.reduce(:+).size
input_size =  tigrNamesInput.values.flatten.uniq.size
total_size =  roleName2TigrNames.values.flatten.uniq.size


##########################################################
if is_output_role
  tigrNamesInput.values.flatten.uniq.each do |tigrName|
    puts [tigrName, tigrName2RoleNames[tigrName]].flatten.join("\t")
  end
  exit
end


##########################################################
results = Parallel.map(tigrNamesInput, in_processes: cpu) do |roleName, tigrNames|
  #puts [roleName, tigrNames.size, roleName2TigrNames[roleName].size].join("\t")
  #puts tigrNames.map{|i|roleName}.join("\n");next
  nums = [tigrNames.size, roleName2TigrNames[roleName].size, input_size, total_size]
  num_str = nums.join(',')
  pvalue = `fisher_chi_test.py --num #{num_str} --type fisher`.chomp
  pvalue_info[roleName] = pvalue.to_f
  pvalue_info
end

results.each do |i|
  pvalue_info.merge!(i)
end

pvalue_info = pvalue_info.sort.to_h


pvalues = pvalue_info.values


fdrs = `fdr_correction0.py  --num #{pvalues.join(',')}`.chomp.split(',').map{|i|i.to_f}
pvalue_info.keys.zip(fdrs) do |roleName, fdr|
  fdr_info[roleName] = fdr
end

fdr_info.each_pair do |roleName, fdr|
  puts [roleName, pvalue_info[roleName], fdr].join("\t")
end


