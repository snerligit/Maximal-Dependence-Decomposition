#! /usr/perl/bin -w

open(READ_FILE_1, "<",$ARGV[0]) or die("Couldn't open file:$!");
open(READ_FILE_2, "<",$ARGV[1]) or die("Couldn't open file:$!");
my @arr;
while(my $line = <READ_FILE_1>)
{
	chomp($line); 
	push(@arr, split("\t",$line));
}
if(scalar(@arr) == 36)
{
	$limit = 9;
}
else
{
	$limit = 14;
}
my @arr_ACGT = @arr;
splice(@arr);

open(WRITE_FILE, ">output.out") or die("Couldn't open file:$!");

while($line = <READ_FILE_2>)
{
	# extract introns or exons depending on the threshold.
	chomp($line);
	push(@arr, split('',$line));
}

$threshold = 0;
# logic to determine the exons and introns and print them out.
for(my $i=0;$i<@arr-$limit+1;$i++)
{
	my @new_arr = @arr;
	my @sub_seq = splice(@new_arr,$i,$limit);
	# for 5' splice site
	if( $limit == 9 && $arr[$i+3] eq 'G' && $arr[$i+4] eq 'T')
	{
		my $score = calc_score($limit,@arr_ACGT,@sub_seq);
		if($score > $threshold)
		{
			printf(WRITE_FILE "$i\t@sub_seq\t%.3f\n",$score);
		}
	}
	# for 3' splice site
	if( $limit == 14 && $arr[$i+10] eq 'A' && $arr[$i+11] eq 'G')
	{
		my $score = calc_score($limit,@arr_ACGT,@sub_seq);
		if($score > $threshold)
		{
			printf(WRITE_FILE "$i\t@sub_seq\t%.3f\n",$score);
		}
	}
}

close(READ_FILE_1) or die("Couldn't close file:$!");
close(READ_FILE_2) or die("Couldn't close file:$!");
close(WRITE_FILE) or die("Couldn't close file:$!");

sub calc_score
{
	my $start = 0;
	my $limit = shift(@_);
	my @arr_A = splice(@_,$start,$limit);
	my @arr_C = splice(@_,$start,$limit);
	my @arr_G = splice(@_,$start,$limit);
	my @arr_T = splice(@_,$start,$limit);
	my @seq = splice(@_,$start,$limit);
	
	my $sum = 0;
	for(my $i=0;$i<@seq;$i++)
	{
		if($seq[$i] eq 'A')
		{
			$sum = $sum + $arr_A[$i];
		}
		else
		{
			if($seq[$i] eq 'C')
			{
				$sum = $sum + $arr_C[$i];
			}
			else
			{
				if($seq[$i] eq 'G')
				{
					$sum = $sum + $arr_G[$i];
				}
				else
				{
					if($seq[$i] eq 'T')
					{
						$sum = $sum + $arr_T[$i];
					}
				}
			}
		}
	}
	return $sum;
}
exit;