#! /usr/bin/env ruby


#############################################################
require 'getoptlong'
require 'parallel'

require 'Dir'


#############################################################
infile = nil
cpu = 1
outdir = nil
is_force = false

out_fhs = Hash.new


#############################################################
def getProtAccn(infile)
  #<h4 id="mrnaandproteins-1834390">mRNA and Protein(s) </h4>
  #<ol>
  #  <li>
  #    <p><a href="/protein/YP_002826335.1">YP_002826335.1</a> cob(I)yrinic acid a,c-diamide adenosyltransferase [Sinorhizobium fredii NGR234]</p>
  prot_accn = nil
  in_fh = File.open(infile, 'r')
  lines = in_fh.readlines.map{|i|i.chomp}

  is_discontinued = lines.count{|i| i =~ /discontinued/} >= 1 ? true : false

  lines.each_with_index do |line, index|
    line.chomp!
    if is_discontinued
      #<li>protein: <a href="/protein/WP_012556018.1/">WP_012556018.1</a></li>
      if line =~ /protein: \<a href=\"\/protein\/([^\/]+)\/\"/
        prot_accn = $1
        puts prot_accn
        break
      end
    else
      if line =~ /mRNA and Protein\(s\)/
        line = lines[index+3]
        if line =~ /\/protein\/ ([^\\]+) \"/x
          prot_accn = $1
        end
        break
      end
    end
  end

  in_fh.close
  return(prot_accn)
end


#############################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--force', GetoptLong::NO_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
    when '--cpu'
      cpu = value.to_i
    when '--outdir'
      outdir = value
    when '--force'
      is_force = true
  end
end


#############################################################
mkdir_with_force(outdir, is_force)
gene_outdir = File.join(outdir, 'gene')
seq_outdir = File.join(outdir, 'seq')
mkdir_with_force(gene_outdir, is_force)
mkdir_with_force(seq_outdir, is_force)

ok_list_file = File.join(outdir, 'pnas.ok.list')
bad_list_file = File.join(outdir, 'pnas.bad.list')


#############################################################
in_fh = File.open(infile, 'r')
lines = in_fh.readlines.map{|i|i.chomp}
in_fh.close

out_fhs[:ok] = File.open(ok_list_file, 'w')
out_fhs[:bad] = File.open(bad_list_file, 'w')


#############################################################
Parallel.map(lines, in_process:cpu) do |line|
  line_arr = line.split("\t")
  gene, id = line_arr[0, 2]
  if gene =~ /( [^()]+ ) \(/x
    gene = $1
  end

  output = `esearch -db gene -query #{id} 2>/dev/null | efetch -format tabular 2>/dev/null | grep -i #{id}`
  arr = output.chomp.split("\n").map{|i|i.chomp}
  n = arr.size
  if n == 0
    out_fhs[:bad].puts [gene, id, n].join("\t")
    next
  elsif n >= 2
    arr = arr.reject{|i| i.split("\t")[6] =~ /,/}
    n = arr.size
    if n != 1
      out_fhs[:bad].puts [gene, id, n].join("\t")
      next
    end
  end

  out_fhs[:ok].puts [gene, id, n].join("\t")

  accn = arr[0].split("\t")[2]
  gene_outfile = File.join(gene_outdir, id + '.html')
  `wget https://www.ncbi.nlm.nih.gov/gene/#{accn} -q -O #{gene_outfile}`

  prot_accn = getProtAccn(gene_outfile)

  seq_outfile = File.join(seq_outdir, id+'.protein')

  `esearch -db protein -query "#{prot_accn}[accn]" | efetch -format fasta > #{seq_outfile}`
end


#############################################################
out_fhs.each_pair do |type, out_fh|
  out_fh.close
end


