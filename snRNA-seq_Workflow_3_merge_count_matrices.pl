#!usr/bin/perl
#### Usage: perl merged_count_matrix_all_cells_individuals.pl Gene_list(one gene each row) output_prefix ####################
my $out = $ARGV[1];


open IN,"< $ARGV[0]" || die "$!";

my %ga;
my @ga;
while(<IN>)
{
    chomp;
    my @tmp = split /\t/,$_;
    if($tmp[0]=~/\bx\b/)
    {}
    else{
        $ga{$tmp[0]} = $tmp[0];
        push @ga,$tmp[0];
    }
}

my @file = @ARGV = <exMatrix_*>;
my $all = scalar(@file);
my $i = 0;
my %ge;
foreach my $file(@file)
{
    if($i == 0)
    {
        $i++;
        open IN1,"< $file" || die "$!";
        open OUT,">> $out.$i" || die "$!";
        my $len = 0;
        while(<IN1>)
        {
            chomp;
            if(/batch\d+?/)
            {
                print OUT "Gene_id\t$_\n";
            }
            else
            {
                my @tmp = split /\t/,$_;
                if(scalar(@tmp) - 1 > $len)
                {
                    $len = scalar(@tmp) - 1;
                }
                my $pri = $tmp[1];
                foreach(2..(scalar(@tmp)-1))
                {
                    $pri .= "\t$tmp[$_]";
                }
                $ge{$tmp[0]} = $pri;
            }
        }
        foreach(keys %ga)
        {
            if($ge{$_})
            {
                $ga{$_} .= "\t$ge{$_}";
            }
            else
            {
                foreach my $num(1..$len){
                    $ga{$_} .= "\t0";
                }
            }
        }
        foreach(@ga)
        {
            print OUT "$ga{$_}\n";
        }
        %ge = undef;
        undef %ge;
    }
    else
    {
        my $out1 = "$out.$i";#file output last time
        my $j = $i + 1;
        my $out2 = "$out.$j";
        print "$out1\t$out2\n";
        open IN1,"< $file" || die "$!";#read the file from last output
        open IN2,"< $out1" || die "$!";#read the file we want to add with
        open OUT,">> $out2" || die "$!";#read file output update

        my $head;
        my $len = 0;
        while(<IN1>)
        {
            chomp;
            my @tmp = split /\t/,$_;
            $len = scalar(@tmp) - 1 unless $len > scalar(@tmp) - 1;
            if(/batch\d+?/)
            {
                $head = $_;
            }
            else
            {
                my $pri = $tmp[1];
                    foreach(2..(scalar(@tmp)-1))
                    {
                        $pri .= "\t$tmp[$_]";
                    }
                    $ge{$tmp[0]} = $pri;
            }
        }
        while(<IN2>){
            chomp;
            my @tmp = split /\t/,$_;
            if(/Gene_id/){
                print OUT "$_\t$head\n";
            }
            else
            {
                if($ge{$tmp[0]})
                {
                    print OUT "$_\t$ge{$tmp[0]}\n";
                }
                else
                {
                    print OUT "$_";
                    foreach(1..$len)
                    {
                        print OUT "\t0";
                    }
                    print OUT "\n";
                }
            }
        }
        %ge = undef;
        undef %ge;
        unlink "$out1" unless $i >= $all - 2;
        $i++;
    }
}
