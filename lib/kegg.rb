#! /usr/bin/env ruby


def getGene2Ortho(gene2OrthoFile, include_kos=[])
  gene2ortho = Hash.new
  in_fh = File.open(gene2OrthoFile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    gene, ortho = line.split("\t")
    ortho = ortho.split(':')[-1]
    if not include_kos.empty?
      next if not include_kos.include?(ortho)
    end
    gene2ortho[gene] = ortho
  end
  in_fh.close
  return(gene2ortho)
end


