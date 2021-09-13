
$input="libs/zm_star0_se.json";

$output = json_read($input);
%output = %{$output};
#print keys %output;
foreach (keys %output) {
	print "$_\t$output{$_}\n";
}

print "\n";

###############################################
# module: preset_read
###############################################
sub json_read {
# read preset file
	my $inpreset_file = shift @_;
	my %outpreset_para;
	my $para_name;
	my $para_content;
	open(PRESET, "<", $inpreset_file) || die;
	while (<PRESET>) {
		chomp;
		if (/^\{$/) {
			do {
				$_ = <PRESET>; chomp;# go next
				&skip_comment;
				if (/\"(\S+)\" *: *\"(.*)\"/) {
					$para_name = $1;
					$para_content = $2;
				} elsif (/\"(\S+)\" *: *\{/) {
					$para_name = $1;
					my @para_content = ();
					do {
						$_ = <PRESET>; chomp; # go next
						&skip_comment;
						if (/\"(\S+)\" *: *\"(.+)\"/) {
							my $subpara = $1." ".$2;
							push(@para_content, $subpara);			
						}
					} until (/^\}$/);
					$para_content = join(" ", @para_content);
				}

				$outpreset_para{$para_name} = $para_content;
			} until (/^\}$/);
		}
	}
	close PRESET;
	if (%outpreset_para) {
		return \%outpreset_para;
	} else {
		#&message_print_quit("no parameters in $inpreset_file"); 
	}

	###############################
	sub skip_comment {
	# module to skip comments
		if (/\"_comment\"/) {
			$_ = <PRESET>; chomp;
		}
	}
	###############################
}


