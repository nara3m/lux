open A, "$ARGV[0]";
my $i=0;
$sum=0;

while (<A>){
chomp ($_);
@line = split (" ", $_);

$line[0] =~ s/H//ig;			# Hydrogen atoms are removed
$line[0] =~ tr/[+,\-,*,.]//d;		# Charges (+,-), wild character (*) and join character (.) are renoved
$line[0] =~ s/\[Na\]//g;			# Sodium ions are removed
$line[0] =~ tr/1-9//d;			# Numbers 1-9 are removed
$line[0] =~ tr/[]//d;			# Square brackets are removed
$line[0] =~ s/Br/A/g;			# Two character atoms (Br, Cl) are given single character
$line[0] =~ s/Cl/D/g;
$line[0] =~ tr/()//;		# Round brackets (meaning - branch opening and closing) are given K and L characters
$line[0] =~ s/\\//g;		# Forward an dbackward slashes (meaning - cis - trans) are given M and Q characters
$line[0] =~ s/\///g;

#print "$line[0]\n";print length($line[0]);<stdin>;
$sum=$sum+length($line[0]);
$i++;
}

close A;

$avg_length = int($sum/$i);

print "i=$i\navg length=$avg_length\n";

__END__
