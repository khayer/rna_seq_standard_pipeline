# Convert bedpe file to interact format (see ucsc http://genome.ucsc.edu/goldenPath/help/interact.html)
# by Katharina Hayer 8/2/21 
# Usage: ruby convert_bedpe_to_interact.rb in.bedpe out.interact skip kind 
# skip == how many lines to skip in in.bedpe (omit header)
# kind == either mustache or cLoops

# e.g
# cd /Users/hayerk/Dropbox_Work/files/BassingLab/2021_11_Bassing_Oltz_Writing_Weekend/merge_HiCCUPS
# ruby ~/github/hic_pipe/scripts/convert_bedpe_to_interact.rb DP_all_reps_hiccups/merged_loops.bedpe DP_all_reps_hiccups.interact 2 hiccups
puts ARGV
outfile = File.open(ARGV[1], 'w')
skip = ARGV[2].to_i
kind = ARGV[3] 
name = "kat_fun"
color = 0 # Use 0 and spectrum setting to shade by score


if kind == "cLoops"
	File.open(ARGV[0]).each do |line|
		line.chomp!
		if (line =~ /^#/) or skip != 0
			skip = skip - 1
			next
	
		end
		fields = line.split("\t")
		out = "#{fields[0]}\t#{fields[1]}\t#{fields[5]}\t#{fields[8]}\t"
		score_double = fields[10].to_f # Using the enrichment score
		score_int = -1
		case
		when score_double > 5
			score_int = 955
		when score_double > 3
			score_int = 722
		when score_double > 2.5
			score_int = 388
		else 
			score_int = 0
		end
	
		
		out += "#{score_int}\t#{score_double}\t#{name}\t#{color}\t"
		source = "#{fields[0]}\t#{fields[1]}\t#{fields[2]}\tAnchorA\t."
		target = "#{fields[3]}\t#{fields[4]}\t#{fields[5]}\tAnchorB\t."
		out += "#{source}\t#{target}"
		#puts line
		outfile.puts out
		#break
	end
elsif kind == "hiccups"
	File.open(ARGV[0]).each do |line|
		line.chomp!
		if (line =~ /^#/) or skip != 0
			skip = skip - 1
			next
	
		end
		fields = line.split("\t")
		out = "#{fields[0]}\t#{fields[1]}\t#{fields[5]}\t#{fields[8]}\t"
		score_double = (-Math.log10(fields[17].to_f+4.034E-50)).round(4) #fields[17].to_f # Donut
		score_int = -1
		case
		when score_double > 5
			score_int = 955
		when score_double > 3
			score_int = 722
		when score_double > 2.5
			score_int = 388
		else 
			score_int = 0
		end
	
		
		out += "#{score_int}\t#{score_double}\t#{name}\t#{color}\t"
		source = "#{fields[0]}\t#{fields[1]}\t#{fields[2]}\tAnchorA\t."
		target = "#{fields[3]}\t#{fields[4]}\t#{fields[5]}\tAnchorB\t."
		out += "#{source}\t#{target}"
		#puts line
		outfile.puts out
		#break
	end
	
else
	count = 0
	File.open(ARGV[0]).each do |line|
		line.chomp!
		if (line =~ /^#/) or skip != 0
			skip = skip - 1
			next
	
		end
		fields = line.split("\t")
		out = "chr#{fields[0]}\t#{fields[1]}\t#{fields[5]}\tmustache#{count}\t"
		score_double = (-Math.log10(fields[6].to_f)).round(4) # Using the FDR
		score_int = -1
		case
		when score_double > 2
			score_int = 955
		when score_double > 1.5
			score_int = 722
		when score_double > 1
			score_int = 388
		else 
			score_int = 0
		end
	
		
		out += "#{score_int}\t#{score_double}\t#{name}\t#{color}\t"
		source = "chr#{fields[0]}\t#{fields[1]}\t#{fields[2]}\tAnchorA\t."
		target = "chr#{fields[3]}\t#{fields[4]}\t#{fields[5]}\tAnchorB\t."
		out += "#{source}\t#{target}"
		#puts line
		outfile.puts out
		count += 1
		#break
	end

end


#DESIRED output

#table interact
#"interaction between two regions"
#    ( 
#    string chrom;        "Chromosome (or contig, scaffold, etc.). For interchromosomal, use 2 records" 
#    uint chromStart;     "Start position of lower region. For interchromosomal, set to chromStart of this region" 
#    uint chromEnd;       "End position of upper region. For interchromosomal, set to chromEnd of this region"
#    string name;         "Name of item, for display.  Usually 'sourceName/targetName/exp' or empty"
#    uint score;          "Score (0-1000)"
#    double value;        "Strength of interaction or other data value. Typically basis for score"
#    string exp;          "Experiment name (metadata for filtering). Use . if not applicable"
#    string color;        "Item color.  Specified as r,g,b or hexadecimal #RRGGBB or html color name, as in //www.w3.org/TR/css3-color/#html4. Use 0 and spectrum setting to shade by score"
#    string sourceChrom;  "Chromosome of source region (directional) or lower region. For non-directional interchromosomal, chrom of this region."
#    uint sourceStart;    "Start position in chromosome of source/lower/this region"
#    uint sourceEnd;      "End position in chromosome of source/lower/this region"
#    string sourceName;   "Identifier of source/lower/this region"
#    string sourceStrand; "Orientation of source/lower/this region: + or -.  Use . if not applicable"
#    string targetChrom;  "Chromosome of target region (directional) or upper region. For non-directional interchromosomal, chrom of other region"
#    uint targetStart;    "Start position in chromosome of target/upper/this region"
#    uint targetEnd;      "End position in chromosome of target/upper/this region"
#    string targetName;   "Identifier of target/upper/this region"
#    string targetStrand; "Orientation of target/upper/this region: + or -.  Use . if not applicable"
#
#    )