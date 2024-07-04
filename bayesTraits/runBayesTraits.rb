#! /usr/bin/env ruby


##########################################################
dir = File.dirname($0)
$: << File.join(dir, "../phyloGLM/")


##########################################################
require 'getoptlong'
require 'parallel'

require 'Dir'
require 'runPhyloglm'


##########################################################
$is_cover = true
$is_varRates = false

indir = nil
treefile = nil
infiles = Array.new
outfile = nil
outdir = nil
cpu = 1
range = Array.new
min_count_gene1 = 1
is_force = false
is_tolerate = false
types = Array.new
stones = [100, 1000]
restricts = Array.new
is_bayestraits = true
bayestraits_args = Hash.new


phyloglm_with_type_info = Hash.new{|h,k|h[k]={}}


##########################################################
class Phyloglm
  attr_accessor :lr_statistic, :pvalue, :cov, :fdr, :bonferroni, :transition_rate, :transition_rate_infos, :restrictPvalueInfo, :restrictQvalueInfo, :restrictBonferroniInfo
  def initialize(gene)
    @gene = gene
    @pvalue = nil
    @transition_rate = Hash.new
    @restrictPvalueInfo = Hash.new
    @restrictQvalueInfo = Hash.new
    @restrictBonferroniInfo = Hash.new
  end
end


##########################################################
def check_requirements()
  %w[chi2 BayesTraitsV3].each do |prog|
    `which #{prog} 2>/dev/null`
    if $?.exitstatus != 0
      raise "Fatal error!\n#{prog} was not installed! Exiting ......"
    end
  end
end


def selectInfilesBasedOnRange(infiles, range)
  if ! range.empty?
    infiles = infiles[range[0]-1, range[1]-range[0]+1]
  end 
  return(infiles)
end


def filterUniValue(infiles, min_count_gene1=1)
  new_infiles = Array.new
  infiles.each do |infile|
    items = Array.new
    count_gene1 = 0
    in_fh = File.open(infile, 'r')
    in_fh.each_line do |line|
      line.chomp!
      line_arr = line.split("\t")
      items << line_arr[2]
      count_gene1 += 1 if line_arr[2] == '1'
    end
    in_fh.close

    if items.uniq.size >= 2 and count_gene1 >= min_count_gene1
      new_infiles << infile
    end
  end
  return(new_infiles)
end


def outputCmdfile(outdir, type, corename, stones=[], add_args={})
  cmdfiles = Array.new
  dep_indep_nums = %w[2 3]
  outdir1 = File.join(outdir, type, corename) # e.g., outdir/MCMC/domain1
  case type
    when /ml/i
      method_num = '1'
    when /mcmc/i
      method_num = '2'
  end

  dep_indep_nums.each do |dep_indep_num|
    outdir2 = File.join(outdir1, dep_indep_num+'-'+method_num) # e.g., outdir/MCMC/domain1/2-3
    mkdir_with_force(outdir2, false, true)
    cmdfile = File.join(outdir2, 'cmdfile')
    cmdfiles << cmdfile
    out_fh = File.open(cmdfile, 'w')
    out_fh.puts dep_indep_num
    out_fh.puts method_num
    out_fh.puts "setMinMaxRate 0 #{add_args[:max_rate]}" if add_args.include?(:max_rate)
    if type =~ /mcmc/i
      out_fh.puts 'varRates' if $is_varRates
      #out_fh.puts "Iterations 100000" #Iterations
      out_fh.puts ["Stones", stones].flatten.join(' ')
    end
    out_fh.puts "LF #{outdir2}/#{corename}"
    out_fh.puts 'run'
    out_fh.close
  end

  return(cmdfiles)
end


def outputCmdfile_restrict(outdir, type, corename, restricts=[])
  cmdfiles = Array.new
  dep_indep_nums = %w[3]
  outdir1 = File.join(outdir, type, corename) # e.g., outdir/MCMC/domain1
  method_num = '1'

  dep_indep_nums.each do |dep_indep_num|
    restricts.each do |restrict| #restrict is an array, e.g. %w[q12 q34 q13]
      outdir2 = File.join(outdir1, dep_indep_num+'-'+method_num, restrict.join('-')) # e.g., outdir/MCMC/domain1/2-3/q12-q34
      mkdir_with_force(outdir2, false, true)
      cmdfile = File.join(outdir2, 'cmdfile')
      cmdfiles << cmdfile
      out_fh = File.open(cmdfile, 'w')
      out_fh.puts dep_indep_num
      out_fh.puts method_num
      restrict_str = restrict.join(' ')
      out_fh.puts "restrict #{restrict_str}"
      #out_fh.puts "restrict q34 q12"
      #out_fh.puts "restrict q24 q13"
      out_fh.puts "LF #{outdir2}/#{corename}"
      out_fh.puts 'run'
      out_fh.close
    end
  end

  return(cmdfiles)
end


def getResultFile(cmdfile, type)
  dir = File.dirname(cmdfile)
  infile = ''
  case type
    when /ml/i
      infile = Dir.glob("#{dir}/*Log.txt").shift
    when /mcmc/i
      infile = Dir.glob("#{dir}/*Stones.txt").shift
  end
  return(infile)
end


def parseBayestraitsRes(cmdfiles, type, restricts)
  log_values = Array.new
  transition_rate_infos = Array.new

  cmdfiles.each do |cmdfile|
    transition_rate_info = Hash.new
    infile = getResultFile(cmdfile, type)

    in_fh = File.open(infile, 'r')
    all_lines = in_fh.readlines
    line = all_lines[-1]
    line_arr = line.split("\t")
    log_value = line_arr[1].to_f
    log_values << log_value

    #if not restricts.empty?
    q_names = all_lines[-2].split("\t").select{|i|i=~/^q..$/}
    q_names.each_with_index do |q_name, index| # %w[q12 q13 q21 q24 q31 q34 q42 q43]
      transition_rate_info[q_name] = line_arr[2+index].to_f
    end
    transition_rate_infos << transition_rate_info
    
    #Tree No	Lh	q12	q13	q21	q24	q31	q34	q42	q43	Root - P(0,0)	Root - P(0,1)	Root - P(1,0)	Root - P(1,1)	
    #1	-172.731580	5.354393	2.121383	28.454021	29.617564	10.123757	19.982465	100.000000	7.279124	0.252442	0.250509	0.247251	0.249799	
    #1	-169.188558	5.089900	34.283743	34.283748	5.089901	0.250000	0.250157	0.249843	0.250000	
    #Log marginal likelihood:	-218.826878
    in_fh.close
  end
  return([log_values, transition_rate_infos])
end


def runBayestraits(treefile, infiles, outdir, types, stones, restricts, is_bayestraits, cpu, bayestraits_args)
  phyloglm_with_type_info = Hash.new{|h,k|h[k]={}}
  log_value_info = Hash.new{|h,k|h[k]={}}

  results = Parallel.map(infiles, in_processes:cpu) do |infile|
    phyloglm_with_type_info = Hash.new # This phyloglm_with_type_info is different from the one outside the Parrallel loop.
    corename = File.basename(infile)
    types.each do |type|
      phyloglm_with_type_info[type] = Hash.new
      cmdfiles = outputCmdfile(outdir, type, corename, stones, bayestraits_args)
      performBayestraits(cmdfiles, treefile, infile, type) if is_bayestraits

      log_values, transition_rate_infos = parseBayestraitsRes(cmdfiles, type, restricts)
      #log_value_info[type][corename] = 2 * log_values.reduce(:-).abs
      log_value_info[type][corename] = 2 * log_values.reverse.reduce(:-)
      #puts [infile, log_value_info[type][corename]].join("\t")
      pvalue = getPvalue(type, log_value_info, corename)

      phyloglm = Phyloglm.new(corename)
      phyloglm.lr_statistic = log_value_info[type][corename]
      phyloglm.pvalue = pvalue
      phyloglm.transition_rate_infos = transition_rate_infos

      phyloglm.restrictPvalueInfo = getRestrictInfo(restricts, is_bayestraits, outdir, type, corename, treefile, infile, log_values)
      phyloglm.restrictPvalueInfo.each_pair do |res, pvalue|
        phyloglm.restrictBonferroniInfo[res] = [pvalue * infiles.size, 1].min
      end

      transition_rate_info = phyloglm.transition_rate_infos[1]
      phyloglm.transition_rate['all'] = %w[q21 q34 q24 q31].map{|i|transition_rate_info[i]}.reduce(:+) - %w[q12 q43 q42 q13].map{|i|transition_rate_info[i]}.reduce(:+)
      if not restricts.empty?
        phyloglm.transition_rate['34-12'] = %w[q34].map{|i|transition_rate_info[i]}.reduce(:+) - %w[q12].map{|i|transition_rate_info[i]}.reduce(:+)
        phyloglm.transition_rate['24-13'] = %w[q24].map{|i|transition_rate_info[i]}.reduce(:+) - %w[q13].map{|i|transition_rate_info[i]}.reduce(:+)
      end

      phyloglm_with_type_info[type][corename] = phyloglm
    end
    phyloglm_with_type_info
  end

  results.each do |h1|
    h1.each do |type, h2| 
      phyloglm_with_type_info[type].merge!(h2)
    end
  end

  return(phyloglm_with_type_info)
end


def getRestrictInfo(restricts, is_bayestraits, outdir, type, corename, treefile, infile, log_values0)
  return {} if restricts.empty?
  restrictPvalueInfo = Hash.new
  log_value_info = Hash.new{|h,k|h[k]={}}

  cmdfiles = outputCmdfile_restrict(outdir, type, corename, restricts)
  performBayestraits(cmdfiles, treefile, infile, type) if is_bayestraits

  log_values, transition_rate_infos = parseBayestraitsRes(cmdfiles, type, restricts)
  log_values.each_with_index do |log_value, index|
    restrict = restricts[index]
    log_value_info[type][restrict] = 2 * (log_value-log_values0[1]).abs
  end
  restricts.each do |restrict|
    pvalue = getPvalue(type, log_value_info, restrict, 1)
    restrictPvalueInfo[restrict] = pvalue
  end

  return(restrictPvalueInfo)
end


def performBayestraits(cmdfiles, treefile, infile, type)
  cmdfiles.each do |cmdfile|
    resfile = getResultFile(cmdfile, type)
    if not resfile.nil? and File.exists?(resfile) and not $is_cover
      next
    else
      system("BayesTraitsV3 #{treefile} #{infile} < #{cmdfile} >/dev/null")
    end
  end
end


def getPvalue(type, log_value_info, corename, fd=4)
  if type == 'ml'
    a = `chi2 #{fd} #{log_value_info[type][corename]}`
    if a =~ /invalid/
      pvalue = 1
    else
      a.chomp!.chomp!
      pvalue = a.split(/\s+/)[-1].to_f
    end
  elsif type == 'mcmc'
    pvalue = log_value_info[type][corename]
  end
  return(pvalue)
end


def get_restrict_fdrs(phyloglm_info, restricts)
  restricts.each do |restrict|
    pvalues_str = phyloglm_info.sort.to_h.values.map{|i|i.restrictPvalueInfo[restrict]}.compact.join(',')
    fdrs = `fdr_correction0.py --num #{pvalues_str}`.chomp.split(',').map{|i|i.to_f}
    phyloglm_info.sort.to_h.each_key do |gene|
      phyloglm = phyloglm_info[gene]
      a = phyloglm.restrictPvalueInfo[restrict]
      a = a.to_s =~ /\d/ ? fdrs.shift : 'NaN'
      phyloglm_info[gene].restrictQvalueInfo[restrict] = a
    end
  end
  return(phyloglm_info)
end



##########################################################
##########################################################
if __FILE__ == $0
  opts = GetoptLong.new(
    ['-t', '--tree', GetoptLong::REQUIRED_ARGUMENT],
    ['-i', '--in', GetoptLong::REQUIRED_ARGUMENT],
    ['-o', '--out', GetoptLong::REQUIRED_ARGUMENT],
    ['--indir', GetoptLong::REQUIRED_ARGUMENT],
    ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
    ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
    ['--range', GetoptLong::REQUIRED_ARGUMENT],
    ['--min_count', GetoptLong::REQUIRED_ARGUMENT],
    ['--force', GetoptLong::NO_ARGUMENT],
    ['--tole', '--tolerate', GetoptLong::NO_ARGUMENT],
    ['--type', GetoptLong::REQUIRED_ARGUMENT],
    ['--stones', GetoptLong::REQUIRED_ARGUMENT],
    ['--restrict', GetoptLong::REQUIRED_ARGUMENT],
    ['--max_rate', GetoptLong::REQUIRED_ARGUMENT],
    ['--var_rates', '--varRates', GetoptLong::NO_ARGUMENT],
    ['--no_c', '--no_cover', GetoptLong::NO_ARGUMENT],
    ['--no_b', '--no_bayestraits', GetoptLong::NO_ARGUMENT],
  )


  opts.each do |opt, value|
    case opt
      when /^(-t|--tree)$/
        treefile = value
      when /^(-i|--in)$/
        infiles << value
      when /^(-o|--out)$/
        outfile = value
      when /^--indir$/
        indir = value
        infiles << read_infiles(value)
      when /^--outdir$/
        outdir = value
      when /^--cpu$/
        cpu = value.to_i
      when /^--range$/
        range = value.split(',').map{|i|i.to_i}
      when '--min_count'
        min_count_gene1 = value.to_i
      when /^--force$/
        is_force = true
      when /^--(tole|tolerate)$/
        is_tolerate = true
      when /^--type$/
        types << value.split(',')
      when /^--stones$/
        stones = value.split(',')
      when /^--restrict$/
        restricts = value.split(',').map{|i|i.split('-')}
      when /^--max_rate$/
        bayestraits_args[:max_rate] = value.to_i
      when /^--(var_rates|varRates)$/
        $is_varRates = true
      when /^--(no_c|no_cover)$/
        $is_cover = false
      when /^--(no_b|no_bayestraits)$/
        is_bayestraits = false
    end
  end


  ##########################################################
  if not is_bayestraits
    raise "Wrong params!\nis_bayestraits is off whereas is_force is on and is_tolerate is off. Exiting ......" if not is_tolerate
  end

  check_requirements()

  types.flatten!

  infiles.flatten!
  infiles = selectInfilesBasedOnRange(infiles, range)
  infiles = filterUniValue(infiles, min_count_gene1)
  puts infiles.size

  mkdir_with_force(outdir, is_force, is_tolerate)

  genes = getAllGenes(indir)

  traitInfos = readDataframeFiles(genes, indir)


  ##########################################################
  phyloglm_with_type_info = runBayestraits(treefile, infiles, outdir, types, stones, restricts, is_bayestraits, cpu, bayestraits_args)

  phyloglm_with_type_info.each_pair do |type, phyloglm|
    phyloglm_info = get_fdrs(phyloglm)
    phyloglm_info = get_restrict_fdrs(phyloglm_info, restricts) if not restricts.empty?
    outfile = File.join(outdir, type) + '.res'
    out_fh = File.open(outfile, 'w')
    out_fh.puts %w[Family LR pvalue qvalue bonferroni no_trait0 prop_trait0 no_trait1 prop_trait1 LR_all LR_res1 LR_res2 p_res1 p_res2 q_res1 q_res2 bonferroni_res1 bonferroni_res2].join("\t")
#TIGR00002	0.66	0.9561	1.0	1	617	1.0	276	1.0	0.0	3.013915	1.0	0.8704	1.0	0.9600163258232235	1	1
    phyloglm_info.sort.to_h.each_pair do |gene, phyloglm|
      traits = traitInfos[gene].sort.to_h.keys
      trait_str = traits.select{|trait|trait!="-"}.map{|trait|[traitInfos[gene][trait].geneNums.reduce(:+), (traitInfos[gene][trait].geneNums.reduce(:+)/traitInfos[gene][trait].species.size.to_f).round(2)].join("\t")}.join("\t")

      transition_rate_str = %w[34-12 24-13].map{|i|phyloglm.transition_rate[i]}.join("\t")
      restrict_fdr_str = ''
      if not restricts.empty?
        restrict_pvalue_str = restricts.map{|restrict|phyloglm.restrictPvalueInfo[restrict]}.join("\t")
        restrict_qvalue_str = restricts.map{|restrict|phyloglm.restrictQvalueInfo[restrict]}.join("\t")
        restrict_bonferroni_str = restricts.map{|restrict|phyloglm.restrictBonferroniInfo[restrict]}.join("\t")
      end
      phyloglm.bonferroni = [phyloglm.pvalue * infiles.size, 1.0].min
      out_fh.puts [gene, phyloglm.lr_statistic.round(2), \
        phyloglm.pvalue, phyloglm.fdr, phyloglm.bonferroni, trait_str, \
        phyloglm.transition_rate['all'], \
        transition_rate_str, restrict_pvalue_str, restrict_qvalue_str, restrict_bonferroni_str].join("\t")
    end
    out_fh.close
  end

end


