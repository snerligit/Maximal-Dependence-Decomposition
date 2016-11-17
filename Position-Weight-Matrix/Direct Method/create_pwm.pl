#!/usr/perl/bin -w

no warnings 'uninitialized';

open (READ_FILE_1, "<",$ARGV[0]) or die("Couldn't open file:$!");
$numerator = 0.1;
$denominator = 0.4;
while($line = <READ_FILE_1>)
{
	chomp($line);
	$row_count++;
	push(@bases, $line);
}

#$row_count = 10;
$den = $row_count+$denominator;
for($i=0;$i<9;$i++)
{
	$count_A = 0;
	$count_C = 0;
	$count_G = 0;
	$count_T = 0;
	for($j=0;$j<$row_count;$j++)
	{
		@each_row = split(//,$bases[$j]);
		if($each_row[$i] eq 'A')
		{
			$count_A++;
		}
		else 
		{
			if($each_row[$i] eq 'C')
			{
				$count_C++;
			}
			else 
			{
				if($each_row[$i] eq 'G')
				{
					$count_G++;
				}
				else
				{
					$count_T++;
				}
				
			}
		}
	}
	push(@arr_A, ($count_A+$numerator)/$den);
	push(@arr_C, ($count_C+$numerator)/$den);
	push(@arr_G, ($count_G+$numerator)/$den);
	push(@arr_T, ($count_T+$numerator)/$den);
}

$exp = 0.25;
for($i=0;$i<@arr_A;$i++)
{
	$arr_A[$i] = (log($arr_A[$i]/$exp))/(log 2);
	$arr_C[$i] = (log($arr_C[$i]/$exp))/(log 2);
	$arr_G[$i] = (log($arr_G[$i]/$exp))/(log 2);
	$arr_T[$i] = (log($arr_T[$i]/$exp))/(log 2);
}
$outfile = $ARGV[0];
$outfile =~ s/_.*//;
$outfile .= "_pwm.out";

open(WRITE_FILE,">",$outfile) or die("Couldn't open file: $!");

print_arr(@arr_A);
print_arr(@arr_C);
print_arr(@arr_G);
print_arr(@arr_T);

sub print_arr
{
	foreach $i(@_)
	{
		printf(WRITE_FILE "%+.3f\t",$i);
	}
	print WRITE_FILE "\n";
	return;
}
close(READ_FILE_1) or die("Couldn't close file:$!");
close(WRITE_FILE) or die("Couldn't close file:$!");
exit;