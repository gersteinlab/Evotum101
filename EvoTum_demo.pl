###package Clonal_fitness3;
#!/usr/bin/perl -w
use lib "/gpfs/home/fas/gerstein/ls926/Scripts/R/Statistics-R-0.33/";
use lib "/gpfs/home/fas/gerstein/ls926/Scripts/R/Statistics-R-0.33/lib/Statistics/";
use lib "/gpfs/home/fas/gerstein/ls926/Scripts/R/Statistics-R-0.33/lib/";
use lib "/gpfs/home/fas/gerstein/ls926/Scripts/R/Regexp-Common-0.01/lib/";
use lib "/gpfs/home/fas/gerstein/ls926/Scripts/R/IPC-Run-0.80/lib/";
use Statistics::R;
use Time::localtime;



$first_name=$ARGV[0];  #file name given by user eg: TESTFILES/
$last_name=$ARGV[1];   #if any, choose a mutation to start from. Helpful if analysis was incomplete to restart from a set point. Often useful to consult $resume_vcf file (see below).
if($first_name)
{
print "Arguments given:\nname-of-file:$first_name\nnr-of-starting-mutation:$last_name\n";
	$run=1;
}else
{
print "Format:perl EvoTum_demo.pl name-of-file [nr-of-starting-mutation]\n";
	$run=0;
}

######### MODES
$fit_expo=1;	#Run the main algorithm. 0 for only generating the intermediate vcf files
$bonfe_mode=1;	#statistical correction based on bonferoni
$driver_list_mode=0; # Downstream analysis included driver detection from coordinates. It does not apply here. 
$resume_mode=1;  # When running a list or multiple jobs, resume_mode=1 deters from running the same vcf again. To rerun delete $resume_file or set to 0
$resume_vcf_mode=1;	# This file keeps track of the last analysed mutation in case the algorightm gets interrupted by R. Then 
$sampling_r_mode=-1;	# Leave '-1' for the default approach with direct Λamda calculation when estimating growth r. 
$generation_mode=1;	# Default measure of k, optimizing for local generational time t_g
$general_classic_mode=1;	# Recommended if interested in first or main clonal mutation. No t_g optimization. Vulnerable to initial perturbation (e.g. CNV). Each hitchiker is generational starting from 1.
$reoptimize_r_mode=0;	# Similar to general_classic_mode, no t_g, but growth r is estimated through NLS optimization. 
$push_mode=1;

############### Sub parameters
$mut_pass_window=100;	#number m of hitchhikers; The size of the sliding window. No k and r values will be reported for the first m sample's mutations. Too large will result in reduced k vales and averaged growth. Too small will be unstable to small changes. 
$coverage=1000;   # Resample Frequency to simulate equal or lower than 1000x coverage
$r_sig_level=0.001;	
#$num_of_ram=50;#	number of random mutations sampled. Not in use in the current demo format.   


##################################### Path and Files
chomp (my @first_name_els=split(/\//,$first_name));
$first_name_wo_path=pop(@first_name_els);  #	separating filename from path
$path=join('/',@first_name_els);  #	path given with file
print "path given is:$path\n";

$pathout="$path"."/OUTFILES/"; #	Self-explanatory; directory of outfiles
$resume_file="$pathout"."RESUME_$first_name_wo_path";	#	This file registers already analyzed files to avoid re-analyzing if running multiple jobs (or a list of files). 
$resume_vcf="$pathout"."RESUME_VCF_$first_name_wo_path";	# This file keeps track of analyzed mutations to help restarting the algorithm from where it stopped. 
$pseudo_vcf_dir="$pathout"."VCF_intermed_$coverage/";# created path for intermidiate files	

#make folders
system ("mkdir $pathout");
system ("mkdir $pseudo_vcf_dir");


######## VCF mode
$vcf_mode=1;
#if(($vcf_mode >0)&&($push_mode eq 0))
#{
#$vcf_file="/gpfs/scratch/fas/gerstein/ls926/TCGA/Pancan/7_mut_SlWindow/Supplemental_Dataset_4-GoldSnvList.tsv";
#$vcf_clone_mode=0;
#$chrom_col=0;
#$pos_col=2;
#$vaf_col=30;
#$ref_nt_col=3;
#$alt_nt_col=4;
#$pat_id_def="300deep_150_limit";
#$pos_end_col=2;

#}elsif(($first_name)&&($last_name)&&($vcf_mode >0)&&($push_mode eq 0))
#{
#$vcf_clone_mode=0;
#$vcf_pseudo="$pseudo_vcf_dir$first_name";
#$vcfid=0;
#$chrom_col=1;
#$pos_col=2;
#$vaf_col=5;
#$ref_nt_col=3;
#$alt_nt_col=4;
#$pat_id_def="$first_name";
#print "resume_vcf:$resume_vcf\tfirst_name:$first_name\n";
#}
if(($vcf_mode >0)&&($push_mode > 0))
{
open RES, "<$resume_vcf";
chomp (my @res_lines=<RES>);
close RES;

$last_line=pop(@res_lines);
chomp (my @last_line_els=split(/\t/,$last_line));
$first_name=$last_line_els[0];
$last_name=$last_line_els[1];

print "fstline_name:$listline_name\tirst_name:$first_name\n";
if($listline_name=~m/$first_name/)
{
	$first_name=$listline_name;

}
else
{
	$first_name=$listline_name;
}


$vcf_clone_mode=0;
$vcf_pseudo="$pseudo_vcf_dir$first_name";
#$vcfid=0;
$chrom_col=0;
$pos_col=1;
$vaf_col=7;
$ref_nt_col=3;
$alt_nt_col=4;
$pat_id_def="$first_name";
print "resume_vcf:$resume_vcf\nvcf_pseudo:$vcf_pseudo\nfirst_name:$first_name\n";
}

$listfile_line="$first_name";
chomp (my @listfile_line_els=split(/_/,$listfile_line));
$listline_name=pop(@listfile_line_els);

#####################

$tm=localtime;
my ($hours,$minutes,$seconds,$day,$month,$year)=($tm->hour,$tm->min,$tm->sec,$tm->mday,$tm->mon,$tm->year);
$fulldate="date$hours-$minutes-$seconds"."_$month-$day-$year";
#print "full date is $fulldate\n";
###############
%mut_freqs=();
%driver_muts_res=();
%driver_muts=();

%success_order=();
%success_k=();
%dri_list=();
%drivers=();
@insim_lines=();
%impact_muts=();
%dri_imp=();
@mutations1=();
%linear_cancers=();

if(($vcf_mode >0)&&($run>0))
{
	if ($resume_vcf_mode>0)
	{
		open VCFLOG, "<$resume_vcf";
		chomp (my @vcflog_lines=<VCFLOG>);
		close VCFLOG;
	
		$vlog_line=pop(@vcflog_lines);
		chomp (my @vlog_line_els=split(/\t/,$vlog_line));
		$pat_test_name=$vlog_line_els[0];			
		$start_sk=pop(@vlog_line_els);
		$start_sk=$start_sk+1;

print "resume vcf is on. Starting from:$start_sk:\n";
	}

	%all_vcf_muts=();	
	%fitness_values_early=();
	@window_muts=();
	@mutations1_parcl=();
	open VCFF, "<$vcf_pseudo";
	chomp (my @vcf_lines=<VCFF>);
	close VCFF;	

	if ($start_sk ge @vcf_lines-1)
	{
system ("rm $resume_vcf");
	}
system ("rm $pseudo_vcf_dir$$pat_id_def");

###############################  Parsing VCF file
	for ($vl=1;$vl<@vcf_lines;$vl++)
	{
		$vcf_line=$vcf_lines[$vl];
		chomp (my @vcf_line_els=split(/\t|_/,$vcf_line));
		$pat_id_formut=$pat_id_def;  # CHECK if pat ID is given in VCF or pseudo VCF format
		$pos_vcf=$vcf_line_els[$pos_col];
		$chrom_vcf=$vcf_line_els[$chrom_col];
		$ref_vcf=$vcf_line_els[$ref_nt_col];
		$alt_vcf=$vcf_line_els[$alt_nt_col];
		$vaf_vcf=$vcf_line_els[$vaf_col];	
		chomp (my @vaf_vcf_els=split(/:|=/,$vaf_vcf));
		$vaf_vcf=pop(@vaf_vcf_els);
		$vaf=2*$vaf_vcf;	# DELETE 2x if frequency is given instead of VAF	

		############ Resampling based on coverage
		if ($coverage >0)
		{
			$vaf_pop=int(1000*$vaf);
			@pop_array=();
			for ($i=0;$i<$vaf_pop;$i++)
			{
				push(@pop_array,1);	
			}
			$remaining_pop=1000-$vaf_pop;
			for ($j=0;$j<$remaining_pop;$j++)
			{
				push(@pop_array,0);
			}
			
			$sum=0;
			for ($c=0;$c<$coverage;$c++)
			{
				$ran_num_gen=int(rand(1000));
				$sum+=$pop_array[$ran_num_gen];				
			}
			$vaf=(int(1000*($sum/$coverage))/1000);				
		}	
	
		$full_mut_vcf="$pat_id_formut"."_$chrom_vcf"."_$pos_vcf"."_$ref_vcf"."_$alt_vcf"."_$vaf";		

		$all_vcf_muts{$full_mut_vcf}=$vaf;

		open PVCF, ">>$pseudo_vcf_dir$pat_id_formut";
print PVCF "$full_mut_vcf\n";
		close PVCF;		
#print "full_mut_vcf:$full_mut_vcf\n";	$include_list_mode
	}
	@mutations1_parcl=sort { $all_vcf_muts{$b} <=> $all_vcf_muts{$a} } keys %all_vcf_muts;
#print "mutations1_parcl:@mutations1_parcl\n";
	$mut_par1=$mutations1_parcl[1];
	$vaftest=$all_vcf_muts{$mut_par1};
	@window_muts=sort { $all_vcf_muts{$b} <=> $all_vcf_muts{$a} } keys %all_vcf_muts;
print "Vaf 1 is  $vaftest. CHECK for prevelance\n";
	%fitness_values_early=CLFIT(\@mutations1_parcl,\@window_muts);
}



sub CLFIT
{

print "clonal_fitness2 mode selected\n";
print "full date is $fulldate\n";
	my ($one_ref, $two_ref)=@_;	
	@mutations1=@{ $one_ref};
	@mutations2=@{ $two_ref};

#print "mutations1:@mutations1\nmutations2:@mutations2\n";

	
	%mut_freqs=();
	%driver_muts_res=();
	%driver_muts=();
	%driver_preds=();
	%success_order=();
	%success_k=();
	%dri_list=();
	%dri_pre=();
	%dri_imp=();
	%drivers=();
	%insims=();
	@insim_lines_ordered=();
	%random_ks=();
	%impact_muts=();
	$resume=1;
#	$run=1;
########### END OF Load prevalent mutations
	$random_mut=$mutations1[4];
	chomp (my @random_mut_els=split(/_/,$random_mut));
	$random_id=$random_mut_els[0];
	if ($vcf_mode >0){$random_id=$pat_id_formut;}




	if ($run>0)
	{

############ CHECKING FOR RESUME MODE
		if (($resume_mode>0)&&($run>0))
		{
			{
				local $/ = undef;
		  		open RFILE, "$resume_file";# or die "Couldn't open file: $!";
		  		binmode RFILE;
		  		$string = <RFILE>;
		  		close RFILE;
		  		if($string=~m/$random_id/){$resume=0;}
			}
		}
		if (($resume>0)&&($run>0))
		{

			@insim_lines=@mutations1;
			$insim_lines_nr=@insim_lines;
			$total_mutations=@mutations2;

			if ($bonfe_mode eq 1)
			{
				$r_sig_level=$r_sig_level/($total_mutations-1);
			}
################### FOR EACH SIMULATION
			$pm=0;
			%insims=();
			%insims_pos=();
			foreach $ins_l(@insim_lines)
			{ 
				chomp (my @ins_l_els=split(/_/,$ins_l));
				$ins_l_freq=pop(@ins_l_els);
				$chpo="$ins_l_els[1]"."_"."$ins_l_els[2]";
				$insims{$ins_l}=$ins_l_freq;
				$insims_pos{$chpo}=$ins_l_freq;
			}

			@insim_lines_ordered=sort({$insims{$b} <=>$insims{$a} }keys %insims);	
#print "insim_lines_ordered:@insim_lines_ordered\n";
			$sl_pl=0;
			$dr_freq=1;	
			%random_muts=();
			$insim_lines_ordered_nr=@insim_lines_ordered;

################### CREATE $num_of_ram Random mutations from sample. TEST IF NEEDED
#			for ($nor=0;$nor<$num_of_ram;$nor++)
#			{
#				$random_num=int rand($insim_lines_nr);
#print "random_num:$random_num\tinsim_lines_nr:$insim_lines_nr\n";
#				$random_muts{$random_num}++;
#				$without_rep=$random_muts{$random_num};
#				if (($without_rep >1)&&($insim_lines_nr >0)&&($insim_lines_nr >= $num_of_ram)){$nor=$nor-1};
#			}
#			@random_muts_keys=keys(%random_muts);
#print "random_muts_keys:@random_muts_keys\n";

			for ($sl=0;$sl<@insim_lines_ordered;$sl++)
			{
				$insim_line=$insim_lines_ordered[$sl];
				chomp (my @insim_line_els=split(/_/,$insim_line));	
				$pat_name=$insim_line_els[0];
				$chrom=$insim_line_els[1];
				$position=$insim_line_els[2];
				$freq=$insim_line_els[5];		
				$mut_freqs{$insim_line}=$freq;
				$mut_test_chrpos="$chrom"."_"."$position";

				if($impact_muts{$mut_test_chrpos})
				{
#			$prev_dri_imp=$dri_imp{$sl};
					$dri_imp{$sl}=$mut_test_chrpos;
				}
		
				if($dr_freq eq $freq)
				{
					$sl_pl=$sl;
#print "New position $sl_pl for freq $freq\n";			
				}else
				{
					if ($sl_pl ne 0)
        				{
        					$dri_list{$sl_pl}++;
						$dri_pre{$sl_pl}="$prediction";
        				}
					$sl_pl=0;
				}
			}

			$dri_list_nr=keys(%dri_list);
#			$driver_list_mode=0;
print "\n$dri_list_nr Preselected drivers:@dri_list_nr_keys\n";
#print @dri_list_nr_keys;
#print "nsim_lines_ordered:@insim_lines_ordered\n";
			@muts_sorted_on_freq=@insim_lines_ordered;

			if (($dri_list_nr >0)||($driver_list_mode eq 0))     ####################### Changed to ge 
			{
				@frequencies=();		
				foreach $fumu(@muts_sorted_on_freq)
				{
					$freq_el=$mut_freqs{$fumu};
					push (@frequencies,$freq_el);
#print "$freq_el\n";
				}	

#print "@frequencies\n";
	################### end of FOR EACH SIMULATION.PARSE MUTATIONS 		
				@to_fit_parent_sorted=@frequencies;  ###CHANGE IF UNSORTED!!!!
				$nr_of_sorted_keys=@to_fit_parent_sorted;
				@ks=();
				%order_tests=();
				%significant_drivers=();
				%p_values_drivers=();
				$dr=0;
				$k_pval=1;
				%drivers=();
				$kval_mode=0;
	################### Identify new drivers	###############	 


				%sk_lines=();
				for($sk=$mut_pass_window;$sk<@frequencies;$sk++)
				{
					$run=1;
					open LOGOUT, ">>$resume_vcf";
print LOGOUT "$random_id\t$sk\n";
					close LOGOUT;	
		
					if($sk eq @frequencies-1)
					{
						system ("rm $resume_vcf");
					}

					if ($sk eq $mut_pass_window)
					{
						open ORFILE, ">>$resume_file";
print ORFILE "$random_id\n";
			                	close ORFILE;
					}	
	
					@r2print=();
					@k2print=();	
#		$sk=$sk-1;    ################ DELETE WHEN IMPROVE CANDIDATES!!!!!!!!!!!!!!!!!!
					$k_pval=1;
			                $dr=0;
			                $mut_to_fit=$to_fit_parent_sorted[$sk]; ################ change to $sk  WHEN IMPROVE CANDIDATES!!!!!!!!!!!!!!!!!!
			                @passenger_mutations=();
			                %pm_freq=();
					$full_mut=$insim_lines_ordered[$sk];
#			                for ($pm=0;$pm<$sk;$pm++)
#			                {
#			                       $passenger_mutation=$to_fit_parent_sorted[$pm];
#			                       $passenger_mutation_frequency=$passenger_mutation;
#			                       push (@passenger_mutations,$passenger_mutation_frequency);
#			                       $pm_freq{$passenger_mutation_frequency}=$pm;
#			                }
			
					$fit_clo_nr=@mutations1;
			                $total_muts=@mutations2;
			                $nr_of_mutations=@mutations2;
			                $fit_tot_nr=@mutations2;
###                        $RFM=$fit_clo_nr-$pm;
					$RFM=$mut_to_fit*$total_muts;
		                        $total_muts=$fit_tot_nr;
		                        $nr_of_mutations=$fit_tot_nr;

		                        $score=0;
		                        @erts=();
		                        @freqs=();
		                        @ert_diffs=();
#                        $nr_ofpat_mut=@passenger_mutations;
		                        $av_ert=0;

		                        $Frt1=$RFM/$total_muts;
#                        $test1=$mut_to_fit/$Frt1;
		                        $order_tests{$sk}=$test1;
		                        $test_tm_sum=0;
		                        @totest=();		

					$ratio=($mut_to_fit*$total_muts)/$RFM;
#print "$sk\t$mut_to_fit\t$ratio\n";
					$score=1;
					$r=2;
					@outputs=();
		                	@rs=();
		                	@ks=();
		                	@av_crzs=();
        		        	$tictoc=0;	
					$prev_var_diff=0;
					$incomplete=0;
					$k=1000;
					$k1=1000;	
					$r_count=0;				
					%k_for_r_values=();
					$crz_average=1000;
					$crz=1000;
					$r_pval=1;		
					$regionfile="";
########################### EStimate r (r_est) based on Frequency equations (and get a p value)
###			if (($sk < $start_sk)){$run=0;}else{$run=1;}

					if (($fit_expo>0)&&($sk >= $mut_pass_window)&&($run>0))
					{
						$r_sum=0;
						$r_sum_count=0;
						@r_ests=();
						@k_ests=();
						$k_average=0;
						$r_sum=0;			
			

						@muts_win=();
		                        	@freqs_win=();
                        			@after_freqs_win=();		
                        			@after_muts_win=();
						$driver_status=0;	
						$regionfile="";

	
##				if (exists $dri_list{$sk}){$driver_status="2:";}
						$fitfr=$to_fit_parent_sorted[$sk];

						@sampled_rs=0;
#				@polynom_coef=();
						if ($sampling_r_mode<0)
						{					
							$first_gen_freq=$to_fit_parent_sorted[$sk-$mut_pass_window-1];
							$test_gen_freq=$to_fit_parent_sorted[$sk-$mut_pass_window-1+($mut_pass_window/2)];		
							$last_gen_freq=$to_fit_parent_sorted[$sk-1];				

							push (@freqs_win,$first_gen_freq,$test_gen_freq,$last_gen_freq);
                        		                push (@muts_win,$sk-$mut_pass_window-1,$sk-$mut_pass_window-1+($mut_pass_window/2),$sk-1);
							if(($test_gen_freq ne $first_gen_freq)&&($test_gen_freq ne $last_gen_freq))
                                		        {
								$lamda=($first_gen_freq-$test_gen_freq)/($first_gen_freq-$last_gen_freq);
                	                		        $zero_coef=1-$lamda;
# -b+-SQRT b2 -4ag )/ 2a 	
								$root1=(1+sqrt(1-4*$lamda*(1-$lamda)))/(2*$lamda);
								$root2=(1-sqrt(1-4*$lamda*(1-$lamda)))/(2*$lamda);
#print "root1:$root1\troot2:$root2\n"; 	
								$nroot1=$root1**(2/$mut_pass_window);
								$nroot2=$root2**(2/$mut_pass_window); 
#print "nroot1:$nroot1\tnroot2:$nroot2\n";
								$secnroot1=$nroot1**(1/2);
								$secnroot2=$nroot2**(1/2);
#print "secnroot1:$secnroot1\tsecnroot2:$secnroot2\n";
					
								$rx1=-log($secnroot1);
								$rx2=-log($secnroot2);
#print "rx1:$rx1\trx2:$rx2\n";					
							}	
						}
						@muts_win_full=();
						@freqs_win_full=();
						for ($skw=$sk-$mut_pass_window;$skw<$sk-1;$skw++)	
#				for ($skw=0;$skw<$sk-1;$skw++)   
						{
							@polynom_coef=();
							push(@muts_win_full,$skw);
							$skw_freq=$to_fit_parent_sorted[$skw];
							push (@freqs_win_full,$skw_freq);
							if($sampling_r_mode>0) 
							{
								$lamda=0;
								@polynom_coef=();
								$first_gen_freq=$to_fit_parent_sorted[$sk-$mut_pass_window-1];
								$last_gen_freq=$to_fit_parent_sorted[$sk-1];
								$test_gen_freq=$to_fit_parent_sorted[$skw];
								if(($test_gen_freq ne $first_gen_freq)&&($test_gen_freq ne $last_gen_freq))
								{
									$lamda=($first_gen_freq-$test_gen_freq)/($first_gen_freq-$last_gen_freq);
									$zero_coef=1-$lamda;
									push(@polynom_coef,$zero_coef);
									for ($fc=$sk-$mut_pass_window;$fc<$skw-1;$fc++)
									{
										push(@polynom_coef,0);
									}
									push(@polynom_coef,-1);
									for($lc=$skw;$lc<$sk-2;$lc++)
									{
										push(@polynom_coef,0);
									}
									push(@polynom_coef,$lamda);
									$polysize=@polynom_coef;
#print "polysize:$polysize\tskw:$skw\tsk:$sk\npolynom_coef:@polynom_coef\n";
									@polynom_coef_rev=reverse(@polynom_coef);
#print "polynom_coef_rev:@polynom_coef_rev\n";

my $R = Statistics::R->new();
									$R->set('pcr',[@polynom_coef]);				
									$R->run(q'library(base)');
									$R->run(q'test<-Re(polyroot(pcr))');
									for ($sol=1;$sol<$mut_pass_window;$sol++)
									{
										$R->set('sol',$sol);
										$R->run(q'rsol <-test[sol]');
										$rsol=$R->get('rsol');
										if($rsol>0)
										{	
											$rx=-log($rsol);
										}else{$rx="NA";}
$R->stop();

									}
						
									$freq_win=$to_fit_parent_sorted[$skw];
									push (@freqs_win,$freq_win);
									push (@muts_win,$skw);
								}
							}
		

							if (($rx2 >= 0)&&($rx2 <2))
							{
my $R = Statistics::R->new();
								$R->set('rx2',$rx2);
								$R->set('freqwin',[@freqs_win]);
								$R->set('freqwin_full',[@freqs_win_full]);
								$R->set('mutwin',[@muts_win]);
								$R->set('mutwin_full',[@muts_win_full]);
								$R->set('Tot',$total_muts);
								$R->set('RFM',$RFM);
								$R->set('fit',$fitfr);
#				$R->run(q'mod <- nls(freqwin ~ exp(a + b * mutwin), start = list(a = 0, b = 0))');
								if ($rx2>0)
								{
									@r2print=();
									@k2print=();
#print "rx2 is greater than 0\n";

if ($generation_mode>0)
{	
#########################   Run with pre-estimated r 
###					$R->run(q'mod <- nls(freqwin ~ exp( -rx2 * (b+mutwin))*(1-a)+a, start = list(a = 0.01, b = 1), control=nls.control(maxiter = 10000000, tol = 1e-04, minFactor = 0.000002,printEval = TRUE, warnOnly = TRUE))');
									$R->run(q'mod <- nls(freqwin_full ~ exp( -rx2 * (b+mutwin_full))*(1-a)+a, start = list(a = 0.01, b = 1), control=nls.control(maxiter = 10000000, tol = 1e-04, minFactor = 0.000002,printEval = TRUE, warnOnly = TRUE))');
									$R->run(q'modpar_a <- coef(mod)[1]');
                        				        	$R->run(q'modpar_b <- coef(mod)[2]');
                        				        	$r_mod_a=$R->get('modpar_a');
                        				        	$r_mod_b=$R->get('modpar_b');
#print "r_mod_b:$r_mod_b\tr_mod_a:$r_mod_a\n";
									$generationtime=$r_mod_b;
                        				                if($generationtime >$sk)	#How to treat big generation times. Default:switch to classic mode (t_g=0)
									{
										$R->run(q'mod <- nls(freqwin_full ~ exp( -rx2 * (mutwin_full))*(1-a)+a, start = list(a = 0.01), control=nls.control(maxiter = 10000000, tol = 1e-04, minFactor = 0.000002,printEval = TRUE, warnOnly = TRUE))');
										$r_mod_a=$R->get('modpar_a');
										$r_mod_b=$sk;
										$generationtime=$r_mod_b;
									}
#print "r_mod_b:$r_mod_b\tr_mod_a:$r_mod_a\n";
				                                	$krz=$total_muts*$fitfr-($r_mod_a*$total_muts);
									if ($krz>0){$kopt=log($total_muts*$fitfr)/log($krz);}else{$kopt="NA";}
				                                	if ($kopt<0){$kopt="NA-";}
				
									push (@r2print,$rx2);
									push (@k2print,$kopt);
}

#########################   Run with pre-estimated r not time estimation with double K
if ($general_classic_mode>0)
{		
									$R->run(q'mod <- nls(freqwin_full ~ exp( -rx2 * (mutwin_full))*(1-a)+a, start = list(a = 0.1), control=nls.control(maxiter = 10000000, tol = 1e-04, minFactor = 0.000002,printEval = TRUE, warnOnly = TRUE))');
			
									$R->run(q'modpar_a <- coef(mod)[1]');
        				                        	$r_mod_a=$R->get('modpar_a');
       					                         	$krza=$total_muts*$fitfr-($r_mod_a*$total_muts);
       					                         	if ($krza>0){$kopta=log($total_muts*$fitfr)/log($krza);}else{$kopta="NA";}
       					                         	if ($kopta<0){$kopta="NA-";}

        				                        	push (@r2print,$rx2);
                				                	push (@k2print,$kopta);
}
					
if($reoptimize_r_mode>0)
{
									$R->run(q'mod <- nls(freqwin_full ~ exp( -b*mutwin_full)*(1-a)+a, start = list(a = 0.01, b = 0.001), control=nls.control(maxiter = 10000000, tol = 1e-04, minFactor = 0.000002,printEval = TRUE, warnOnly = TRUE))');
									$R->run(q'modpar_a <- coef(mod)[1]');
									$R->run(q'modpar_b <- coef(mod)[2]');
									$r_mod_a=$R->get('modpar_a');
									$r_mod_b=$R->get('modpar_b');
									$krz=$total_muts*$fitfr-($r_mod_a*$total_muts);
									if ($krz>0){$kopt=log($total_muts*$fitfr)/log($krz);}else{$kopt="NA";}
									if ($kopt<0){$kopt="NA-";}
									push (@r2print,$r_mod_b);
				                                	push (@k2print,$kopt);	
}
								}

								$line_to_push="sk:$sk\tPreselected Driver?:$driver_status\tgenTime=$generationtime\trx2:$rx2\trx1:$rx1\tkopt:@k2print\tloggrowth:$loggrowth\tr:@r2print\tkopt:@k2print\t$full_mut";

								open FILO, ">>$pathout$pat_name"."_cov$coverage";				
print FILO "$line_to_push\n";
print "$line_to_push\n"; ###  DELETED FILE
#								$line_to_push="sk:$sk\tPreselected Driver?:$driver_status\tloggrowth:$loggrowth\tr:$r_mod_b\tkopt:$kopt\t$full_mut\n";
								$sk_lines{$sk}=$line_to_push;
								close FILO;
$R->stop();

							}
							if ($kopt ne "NA")
							{
								$random_ks{$sk}=$kopt;
							}
						}
						if($driver_list_mode eq 0)
						{
#							for ($ske=0;$ske<$sk-1;$ske++)
#	               					{
#								$frm1=$to_fit_parent_sorted[$ske];
#								$frm2=$to_fit_parent_sorted[$ske+1];
#								$frm3=$to_fit_parent_sorted[$ske+2];
#								$frm_fit2=$to_fit_parent_sorted[$sk];
#								if(($frm1 ne $frm2)&&($frm2 ne $frm3))
#								{
#									$r_sum_count++;	
#									$r_part_est=log(($frm1-$frm2)/($frm2-$frm3));
#									$r_sum+=log(($frm1-$frm2)/($frm2-$frm3));
#									push (@r_ests,$r_part_est);
#								}
#print "k_ests:@k_ests\trealK:$fit_fit\n";
#							}
#			$r_est_average=$r_sum/$r_sum_count;
#							$frm_fit=$to_fit_parent_sorted[$sk+1];
#							$frm_fit1=$to_fit_parent_sorted[$sk-1];
#							$frm_fit2=$to_fit_parent_sorted[$sk];
#							if (($frm_fit ne $frm_fit2)&&($frm_fit1 ne $frm_fit2))#&&($frm1 ne $frm2)&&($frm2 ne $frm3))
#							{
#								$frm_fit2_r_est=log(($frm_fit1-$frm_fit2)/($frm_fit2-$frm_fit));	
#								$nr_of_rests=@r_ests;
#				
#my $R = Statistics::R->new();
#								$R->set('vals',[@r_ests]);						
#								$R->set('val',$frm_fit2_r_est);
#								$R->set('ct', $nr_of_rests);
#                       						$R->run(q'vals_mean<-mean(vals)');
#								$R->run(q'vals_median<-median(vals)');
#			                      			$R->run(q'vals_sd<-sd(vals)');
#        			                		$R->run(q'vals_se<-vals_sd/sqrt(ct)');
#                			        		$R->run(q'p_val<-pnorm(val,mean=vals_mean,sd=vals_se)');
#              				         		$r_pval=$R->get('p_val');
#              				         		$r_win_mean=$R->get('vals_mean');
#								$r_win_median=$R->get('vals_median');
#                			        		@r_vals_back=$R->get('vals');
#$R->stop();
#								if ($r_pval<$r_sig_level)
#								{
#print "MUT $sk+1 is canditate\tr_win_mean:$r_win_mean\tr_median:$r_win_median\tpval:$r_pval\nk_ests:@k_ests\trealK:$fit_fit\n";
#									@k_ests=();
#									$k_sum_forav=0;
#									$k_average=0;
##					for ($skk=$sk-$mut_pass_window;$skk<$sk;$skk++)
#									for ($skk=0;$skk<$sk;$skk++)
#                        			        		{
#										$skk1=$skk+1;
#										$r_for_k=$r_ests[$skk];
#										$frm1=$to_fit_parent_sorted[$skk];
#										$fit2=$to_fit_parent_sorted[$sk];
#										if(($r_win_mean ne 0)&&($skk1>0))
#        			                        	       		{
#                        			        	               	        $crz_est=((($frm1*$total_muts)-$RFM)-(exp($r_win_mean*$skk1*(-1))*($total_muts-$RFM)))/(exp($r_win_mean*$skk1*(-1))-1);
#print "#crz_est:$crz_est\n";
#											if ($crz_est>0)
#											{
#                                	        			       		       	 $k_est=log($fit2*$total_muts)/log($crz_est);
#                                	                			       		 push(@k_ests,$k_est);
#												$k_sum_forav+=$k_est;
#											}	
#                                				        	}
#									}
#									$k_ests_nr=@k_ests;
#									if ($k_ests_nr > 0)
#									{
#										$k_average=$k_sum_forav/$k_ests_nr;
#									}else{$k_average="NA";}
#print ">>k_ests:@k_ests\tk_average:$k_average\trealK:$fit_fit\n";	
#print "k_average:$k_average\trealK:$fit_fit\n";						
#								}
#							}	
						}
###########################################################################

#####################################
#						%test_drivers=();

						if ((($incomplete eq 0)&&(($sk ge $mut_pass_window)&&($r_pval<$r_sig_level)))||(exists $dri_list{$sk}))
						{
							$last_clock=$tictoc-2;
	
							$last=pop(@outputs);
							$output=pop(@outputs);
							$hash_key=$muts_sorted_on_freq[$sk];
	
print "\n***** r:$hash_key*****\n$output\n$last\n";				

							$driver_muts_res{$hash_key}=$last;		
						}
					}
				}
			}
			@random_ks_keys=keys(%random_ks);
			$random_ks_keys_nr=@random_ks_keys;
			$average_random_k=0;
			foreach $random_ks_key(@random_ks_keys)
			{
				$random_ks_val=$random_ks{$random_ks_key};
				$average_random_k+=$random_ks_val/$random_ks_keys_nr;
			}		
#print "Average_random_k:$average_random_k \n";	

			return %driver_muts_res;
		}
	}
}
1;
print "Supercalifragilous\n";


#}
















