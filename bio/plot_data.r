#!/usr/bin/env Rscript 

# About:
# Plot data from tab delimited or csv file 
# In particular, the input file can have an arbitrary number of columns

# Example: 
# Suppose we have a three column file such that:
# $ cat example/file_Rplot.txt 
# X	Y1	Y2
# 0	0	0
# 1	10	2
# 2	10	3
# 3	11	8
# 4	12	7
# 5	15	6

# Useage:
# ./plot_data.r -i example/file_Rplot.txt --prefix example/out_Rplot  --title=MyPlot --width 16 --height 10 --label "X,Y1,Y2" --xaxis "X Value" --yaxis="Y Value" --zero 

# Notes:
# This script has some quirks: it expects a header on the input, and thus ignores the first row
# The label field expects a comma-delimited list of names for each column, although that of the first is necessary but ignored
# It also requires some R packages:
# install.packages("optparse")
# install.packages("ggplot2")

# Script Notes:
# if no shebang, run as:
# R --slave --vanilla --args --help < plot_data.r
# R --slave --vanilla --args -i example_plot.csv --prefix SAAn4165  --title="SAAn4165 Genetic Dist vs Position" --width 16 --height 10 --label1 Ingwavuma --label2 Manzanilla --xaxis Pos --yaxis="Genetic Distance" < plot_data.r

library(ggplot2);
library(optparse);
library(reshape);

args = commandArgs(TRUE);

# inputfile=args[1];			# inputfile
# outputprefix=args[2];			# outputprefix
# samplename=args[3];			# name
# s1=args[4];				# sample1 name
# s2=args[5];				# sample2 name

option_list <- list( 
	make_option(c("-i", "--input"), help="input file"), 
	make_option("--csv", action="store_true", default=FALSE, help="input file is comma delimited"),
	make_option("--prefix", default="out", help="output file prefix"), 
	make_option("--sample", default="sample", help="sample name"), 
	make_option("--title", default="Genetic Distance vs Postion", help="plot title"), 
	make_option("--width", type="integer", default=14, help="width (inches)"), 
	make_option("--height", type="integer", default=8, help="length (inches)"), 
	make_option("--v1", type="integer", default=0, help="vert line at this xpos"), 
	make_option("--v2", type="integer", default=0, help="vert line at this xpos"), 
	make_option("--v3", type="integer", default=0, help="vert line at this xpos"), 
	make_option("--label", default="", help="comma delimited list of column names"), 
	make_option("--label1", default="label1", help="label for first column of data"), 
	make_option("--label2", default="label2", help="label for second column of data"), 
	make_option("--xaxis", default="Position", help="label for x axis"), 
	make_option("--yaxis", default="Genetic Distance", help="label for y axis"), 
	make_option("--scatter", action="store_true", default=FALSE, help="make a scatter plot"),
	make_option("--smooth", action="store_true", default=FALSE, help="add smoothed lines"),
	make_option("--zero", action="store_true", default=FALSE, help="x axis starts from zero")
) 

parser <- OptionParser(usage = "For a 3 column tab delimited input:\n./plot_data.r -i example/file_Rplot.txt --prefix example/out_Rplot  --title=MyPlot --width 16 --height 10 --label \"X,Y1,Y2\" --xaxis X_Value --yaxis=Y_Value --zero\n\nNote: This script has some quirks: it expects a header on the input, and thus ignores the first row\nThe label field expects a comma-delimited list of names for each column, although that of the first is necessary but ignored", option_list=option_list) 

arguments <- parse_args(parser, positional_arguments = TRUE) 
opt <- arguments$options 

# make sure input exists
{
	if ( opt$input == "" )
	{
		print("Error: input file not specified");
		# quit
		q();
	}
	else
	{
		inputfile=opt$input;
		print(paste("input file: ",inputfile,sep=""));

	}
}

outputprefix=opt$prefix;
print(paste("prefix: ",outputprefix,sep=""));

samplename=opt$sample;
print(paste("sample: ",samplename,sep=""));

s1=opt$label1;
print(paste("label1: ",s1,sep=""));

s2=opt$label2;
print(paste("label2: ",s2,sep=""));

figwidth=opt$width;
print(paste("width: ",figwidth,sep=""));

figheight=opt$height;
print(paste("height: ",figheight,sep=""));

figtitle=opt$title;
print(paste("title: ",figtitle,sep=""));

figx=opt$xaxis;
print(paste("xax: ",figx,sep=""));

figy=opt$yaxis;
print(paste("yax: ",figy,sep=""));

print(paste("label: ",opt$label,sep=""));

# print empty line
cat("\n");

# initialize dataframe
df=NULL;

# if input non zero
{
	if ( file.exists(inputfile) )
	{
		if ( opt$csv )
		{
			
			df=read.csv(inputfile, header = TRUE);
		}
		else
		{
			df=read.delim(inputfile, header = TRUE);
		}
	}
	else
	{
		print("Error: input file not found");
		# quit
		q();
	}
}

colnum = length(df[1,]);

{
	# check that label size == number of columns
	if ( opt$label != "" && ( length(df[1,]) != length(strsplit(opt$label, ",")[[1]]) ) )
	{
		print("Error: number of labels do not match number of columns in data frame");
		# quit
		q();
	}
}

# rename cols
{
	if ( opt$label == "" )
	{
		# name cols 1 to number of columns
		names(df) <- seq(1,length(df[1,]),1) 
	}
	else
	{		
		names(df) <- strsplit(opt$label, ",")[[1]]
	}
}

# print(df);

# names(df)<- c("pos","a","b")

# do a "melt" by hand
# copy data frames into 2 new data frames
#df1=df
#df2=df

# delete columns:
#df1$b <- NULL
#df2$a <- NULL

# add column "name" to each:
#df1$sample = s1
#df2$sample = s2

# make header the same:
#colnames(df1)<-c("pos","col","sample")
#colnames(df2)<-c("pos","col","sample")

# row bind:
#df3 <- rbind(df1, df2)

# melt the dataframe into 2 col form:
df3 <- melt(df, id=names(df)[1])

# now plot
p <- ggplot(df3);

{
	if (opt$scatter)
	{
		g1 <- geom_point(aes(x=df3[,1], y=df3[,3], color = variable), size = 1);
	}
	else
	{
		g1 <- geom_line(aes(x=df3[,1], y=df3[,3], color = variable), size = 1);
	}
}

# make label
l <-labs(title=figtitle, x=figx, y=figy)

# change colors
#clr <- scale_color_manual(values=c("blue", "red"))
clr <- scale_color_manual(values=rainbow(colnum - 1))

g2 = NULL
extra0=NULL
extra1=NULL
extra2=NULL
extra3=NULL

{
	if ( opt$smooth ) 
	{
		g2 <- geom_smooth(aes(x=df3[,1], y=df3[,3], color = variable), size = .5);
	}
	if ( opt$zero ) 
	{
		extra0 <- coord_cartesian(xlim=c(0, max(df3[,1])))
	}
	if ( opt$v1 > 0 ) 
	{
		extra1 <- geom_vline(xintercept = opt$v1 )
	}
	if ( opt$v2 > 0 ) 
	{
		extra2 <- geom_vline(xintercept = opt$v2 )
	}
	if ( opt$v3 > 0 ) 
	{
		extra3 <- geom_vline(xintercept = opt$v3 )
	}
}

p + g1 + g2 + l + clr + extra0 + extra1 + extra2 + extra3

print(paste("saving file: ", outputprefix, ".pdf", sep=""));
cat("\n");

ggsave(file=paste(outputprefix, ".pdf", sep=""), width=figwidth, height=figheight, units="in")

# hack:
system("rm Rplots.pdf")
