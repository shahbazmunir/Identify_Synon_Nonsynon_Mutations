#!/usr/bin/perl
use strict;
use warnings;

# Custom script
# User input
my ($SNPFile, $ORFfile) = @ARGV;
if(scalar @ARGV != 2){
	print "\nUsage:\n\t./script.pl AB-A_mergedVar.tab refSeq.Annotation.tab \n\n";
	exit;
}

#my $SNPFile    = ;#"AB-A_mergedVar.tab";
#my $ORFfile    = ;#"Ta.seq.uniq.clean.63.HF.Annotation.tab";
my $Outputfile = $SNPFile;	# "AB-A_mergedVar.mutations.tab";
$Outputfile =~ s/tab/mut.tab/;

#---------------Prepare SNP Dictionary{value => key}------------------
my %SNPList;
my @line;

open(my $WRITE, ">$Outputfile") or die "Unable to open a file '$Outputfile' $!";

open(my $DATA, "<$SNPFile") or die "unable to open file '$SNPFile' $!\n";
while(<$DATA>){
        if($_ =~ m/^#/g){
		print $WRITE $_;
		#print $_;
                next;
        }elsif($_ =~ m/</g){
		next;
	}
	chomp $_;
        @line = split("\t", $_);
	my $allSNPs = join("\t", @line[2..(scalar @line -1)]);
        #print "@line\t$allSNPs\n";
        $SNPList{$line[1]} = $allSNPs;
}

#-----------Process ORF and characterize the SNPs--------------------
open(DATA2, "<$ORFfile") or die "Unable to open a file '$ORFfile' $!";
while(<DATA2>){
	# skip header or UTR region
	if($_ =~ m/^#|UTR/){
		next;
	}

	my @ORF = split("\t", $_);
	my $ORFStart = $ORF[3];
	my $ORFEnd = $ORF[4];
	my @DNASeq = split("", $ORF[6]);
	my @SNPs;

	#print "$ORFStart\t$ORFEnd\n";

	my $index = 0;
	for(my $i = $ORFStart; $i<=$ORFEnd; $i++){#search for SNP within ORF
		$index++;
		

		if(exists $SNPList{$i}){
			#print "$i\n";
			if($SNPList{$i} =~ m/\./g){
				print $WRITE "$i\t$SNPList{$i}\n";
				next;
			}


			#print $WRITE "$i\t$SNPList{$i}\n";
			#print $ORF[7];
			 
			#Split all lines SNPs and store it in an array
			@SNPs = split("\t", $SNPList{$i});
			#print "$i\t@SNPs\n";

			my $diff_aa = ".";
			# process all lines SNPs one by one and look for Synonmous and non-synonmous mutations
			for(my $j=0; $j<(scalar @SNPs); $j++){#process 19 lines
	                        my @modifiedSeq = @DNASeq;
				
				#process comma separated SNPs
				my @arr = split(",", $SNPs[$j]);
				for(my $k=0; $k<(scalar @arr); $k++){
					next if($arr[$k] eq ".");
					#Frame shift mutations - Filter indels
					#if($arr[$k] eq "-" || length $arr[$k] > 1 && ((length $arr[$k]) % 3) != 0){
					if(length $arr[$k] > 1 && ((length $arr[$k]) % 3) != 0){
						$SNPs[$j] = "$SNPs[$j];f";
						next;
					}elsif(((length $arr[$k]) % 3) == 0 || $arr[$k] eq "-"){
						next;
					}
					#print "$modifiedSeq[$index]\t";
					$modifiedSeq[$index] = $arr[$k];
					#print "$arr[$k]\n";
					#print @modifiedSeq; print "\n"; print @DNASeq; print "\n";
					#my $modifiedORF = ReadingFrame(join("", @modifiedSeq));
					my $modifiedORF = translateDNAseq(join("", @modifiedSeq));
					#print "$modifiedORF\n$ORF[7]\n";
					chomp $ORF[7];

		# not efficient way - time was too short so just did that
		my @mod_orf = split("", $modifiedORF);
		my @org_orf = split("", $ORF[7]);
		# Nucleotide Triplet dont encode any amino acid, if it contains char 'N'.
		# 
		if(scalar @mod_orf ne scalar @org_orf){
			$diff_aa = "!";
		}else{
			for (my $z=0; $z< (scalar @mod_orf); $z++){
				if($mod_orf[$z] ne $org_orf[$z]){
					$diff_aa = "$mod_orf[$z],$org_orf[$z]";
				}
			}
		}
					if($modifiedORF eq $ORF[7]){
						$SNPs[$j] = "$SNPs[$j];s";
					}else{
						#print "$ORF[7]\n$modifiedORF\n";
					        $SNPs[$j] = ("$SNPs[$j];n");
					}

				}#end of for loop
			
			}
			$SNPList{$i} = join("\t", @SNPs);
			#print "$SNPList{$i}\n";

			#print "$i\t$SNPList{$i}\t$diff_aa\n";
                        print $WRITE "Ta63c\t$i\t$SNPList{$i}\t$diff_aa\n";
			#print "Ta63c\t$i\t$SNPList{$i}\t$diff_aa\n";
		}#end of if condition
	}
}#end of while loop


#--------------------------Subroutines-----------------------
#Generate Reading frames of a DNA sequences
sub ReadingFrame{
        my $DNA = uc($_[0]);
        my $Triplet = "";
        my $ORF = "";
	my $check = "FALSE";
	my $AminoAcid = "";
	my $ORFStart;
	my $ORFEnd;
	my $seq = "";
	my $Output = "";

	#print "GeneAcc\tFrame\tFrom\tTo\tLength\tDNA_Seq\tORF\n";
        my $readingFrame = "";
        for(my $i=0;$i<length($DNA);$i=$i+3){	#Generate Triplets
                $Triplet = substr($DNA, $i ,3);
		$AminoAcid = Translate($Triplet);	#Find the approptiate amino acid
		if($AminoAcid eq 'M'){	#Start of the ORF
			$check = "TRUE";
			$ORFStart = $i;
		}
		if($check eq "TRUE"){
			$ORF .= $AminoAcid;
			$seq .= $Triplet;
		}
		if($AminoAcid eq '*' && $check eq "TRUE"){	#End of the ORF
			$check = "FALSE";
			$ORFEnd = $i;
		}
		if($ORF ne "" && $check eq "FALSE"){	#Output this ORF
			$Output .= $ORF;
			$ORF = "";
			$seq = "";
		}

        }

	return $Output;
}

sub translateDNAseq{
	my $DNA = uc($_[0]);
	my $Triplet = "";
	my $AminoAcid = "";
	my $res = "";
	for(my $i=0;$i<length($DNA);$i=$i+3){   #Generate Triplets
		$Triplet = substr($DNA, $i ,3);
		#print $i.$Triplet;
		$AminoAcid = Translate($Triplet);
		$res .= $AminoAcid;
	}
	#print "$res\n";
	return $res;
}


#Translate nucleotides into amino acids
sub Translate{
        my $Triplet = $_[0];
        my $aa;
        if($Triplet =~ m/GC[ATGCN]/g){  #Ala / A
                return "A";
        }elsif($Triplet =~ m/TG[CT]/g){ #Cys / C
                return "C";
        }elsif($Triplet =~ m/ATG/g){ #Met / M
                return "M";
        }elsif($Triplet =~ m/GA[TC]/g){ #Asp / D
               return "D";
        }elsif($Triplet =~ m/GA[AG]/g){ #Glu / E
                return "E";
        }elsif($Triplet =~ m/TT[TC]/g){ #Phe / F 
                return "F";
        }elsif($Triplet =~ m/GG[ATGCN]/g){ #Gly / G
                return "G";
        }elsif($Triplet =~ m/CA[TC]/g){ #His / H
                return "H";
        }elsif($Triplet =~ m/AT[ATC]/g){ #Ile / I
                return "I";
        }elsif($Triplet =~ m/AA[AG]/g){ #Lys / K
                return "K";
        }elsif($Triplet =~ m/CT[ATGCN]|TT[AG]/g){ #Leu / L 
                return "L";
        }elsif($Triplet =~ m/AA[TC]/g){ #Asn / N
                return "N";
        }elsif($Triplet =~ m/CC[ATGCN]/g){ #Pro / P
                return "P";
        }elsif($Triplet =~ m/CA[AG]/g){ #Gln / Q
                return "Q";
        }elsif($Triplet =~ m/CG[ATGCN]|AG[AG]/g){ # Arg / R
                return "R";
        }elsif($Triplet =~ m/TC[ATGCN]|AG[CT]/g){ #Ser / S
                return "S";
        }elsif($Triplet =~ m/AC[ATGCN]/g){ #Thr / T
                return "T";
        }elsif($Triplet =~ m/GT[AGCTN]/g){ #Val / V
                return "V";
        }elsif($Triplet =~ m/TGG/g){ #Trp / W
                return "W";
        }elsif($Triplet =~ m/TA[TC]/g){ #Tyr / Y
                return "Y";
        }elsif($Triplet =~ m/TA[GA]|TGA/g){ #*
                return "*";
        }

}

