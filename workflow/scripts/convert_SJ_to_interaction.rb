# Convert bedpe file to interact format (see ucsc http://genome.ucsc.edu/goldenPath/help/interact.html)
# by Katharina Hayer 8/2/21 
# Usage: ruby convert_bedpe_to_interact.rb in.bedpe out.interact skip kind 
# skip == how many lines to skip in in.bedpe (omit header)
# kind == either mustache or cLoops

# e.g
# cd /mnt/isilon/weitzman_lab/projects/lab_rna_seq/results/mapped
# ruby ~/data/tools/rna_seq_standard_pipeline/workflow/scripts/convert_SJ_to_interaction.rb SJ.out.tab SJ.interact 5
puts ARGV
outfile = File.open(ARGV[1], 'w')
cutoff = ARGV[2].to_i
#kind = ARGV[3] 
name = "kat_fun"
color = 0 # Use 0 and spectrum setting to shade by score

skip = 0
count = 0
samplename = ARGV[1].sub("_SJ.out.tab","")

outfile.puts "track type=interact name=\"#{samplename}\" description=\"#{samplename} splicing\" interactDirectional=true maxHeightPixels=200:100:50 visibility=full"

File.open(ARGV[0]).each do |line|
	line.chomp!
	if (line =~ /^#/) or skip != 0
		skip = skip - 1
		next

	end

	fields = line.split("\t")
	next unless fields[0] =~ /^chr/
	strand = "."
	if fields[3] == "1"
		strand = "+"
	elsif fields[3] == "2"
		strand = "-"
	end
	score_double = fields[6].to_i
	out = "#{fields[0]}\t#{fields[1].to_i-5}\t#{fields[2].to_i+4}\tj_#{count}_NR#{score_double.to_i}\t"
	
	next if score_double < cutoff
	score_int = -1
	case
	when score_double > 100
		score_int = 955
	when score_double > 50
		score_int = 722
	when score_double > 10
		score_int = 388
	else 
		score_int = 0
	end

	
	out += "#{score_int}\t#{score_double.to_i}\t#{name}\t#{color}\t"
	source = "#{fields[0]}\t#{fields[1].to_i-5}\t#{fields[1]}\tAnchorA_#{score_double.to_i}\t#{strand}"
	target = "#{fields[0]}\t#{fields[2].to_i-1}\t#{fields[2].to_i+4}\tAnchorB\t#{strand}"
	out += "#{source}\t#{target}"
	#puts line
	outfile.puts out
	count += 1
	#break
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