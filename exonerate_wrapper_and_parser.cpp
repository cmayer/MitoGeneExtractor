#include "exonerate_wrapper_and_parser.hpp"
using namespace std;

bool test_exonerate_binary(const faststring &binary_path)
{
   faststring cmd = binary_path + " -v";
   int err = system(cmd.c_str());

   if (global_verbosity>10)
   {
      cout << "Testing existance of exonerate  binary by calling: " << cmd << "\n";
      if (WEXITSTATUS(err) == 1)
         cout << "Binary exists, so we proceed.";
   }
   return WEXITSTATUS(err) == 1;
}


int run_exonerate(const exonerate_parameters &ep)
{


  if (!test_exonerate_binary(ep.exonerate_binary))
  {
     cerr << "The binary " << ep.exonerate_binary << " does not exist or does not work. Exiting.\n";
     exit(-3);
  }

  return run_exonerate(ep.exonerate_binary, ep.protein_reference, ep.fasta_reads, ep.get_exonerate_vulgar_filename(), ep.get_exonerate_log_filename(), ep.genetic_code_number, ep.frameshift_penalty, ep.max_trials, ep.mode);
}

int run_exonerate(const faststring &exonerate_binary, const faststring &protein_reference, const faststring &fasta_reads,
                  const faststring &exonerate_vulgar_filename, const faststring &exonerate_log_filename,
                  char genetic_code_number, int frameshift_penalty, int max_trials, char mode /*=0*/)
{
  faststring code_string;

  if (genetic_code_number == 1) // The Standard Code (transl_table=1)
    code_string = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 2) // The Vertebrate Mitochondrial Code
    code_string = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";
  else if (genetic_code_number == 3)  // The Yeast Mitochondrial Code
    code_string = "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 4)  // The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
    code_string = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 5)  // The Invertebrate Mitochondrial Code
    code_string = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 6)  // The Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6)
    code_string = "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";


  else if (genetic_code_number == 9)  // The Echinoderm and Flatworm Mitochondrial Code
    code_string = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 10) // The Euplotid Nuclear Code (transl_table=10)
    code_string = "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 11) // The Bacterial, Archaeal and Plant Plastid Code (transl_table=11)
    code_string = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 12) // The Alternative Yeast Nuclear Code (transl_table=12)
    code_string = "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 13) // The Ascidian Mitochondrial Code
    code_string = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 14) // The Alternative Flatworm Mitochondrial Code
    code_string = "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 16) // Chlorophycean Mitochondrial Code
    code_string = "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  else
  {
    cerr << "ERROR: Unknown genetic code number passed to function run_exonerate. Code number: " << genetic_code_number << endl;
    good_bye_and_exit(-20);
  }

  //  system "exonerate --geneticcode FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG --frameshift -9 --query ${dir2}/${gene}.fasta
  // -Q protein --target $dir1/${name}.fas -T dna --model protein2dna --showalignment 0 --showvulgar 1 > $dir3/${gene}_${name}_vulgar.txt 2> $dir3/${gene}_${name}_vulgar.log";

  int system_return_value;
  faststring score_threshold;
  if (global_exonerate_score_threshold != UINT_MAX)
    score_threshold = " -s " + faststring(global_exonerate_score_threshold) + " ";
  faststring command =  exonerate_binary + " --geneticcode " + code_string.c_str() + score_threshold.c_str()
  + " --frameshift " + faststring(frameshift_penalty).c_str()
  + " --query " + protein_reference + " -Q protein --target " + fasta_reads
  + " -T dna --model protein2dna --showalignment 0 --showvulgar 1 1> "
  + exonerate_vulgar_filename + " 2> " + exonerate_log_filename;

  // Special debug mode to examine command string:
  if (mode == 'd')
  {
    cout << "DEBUG: Exonerate command: " << command << endl;
    return 0;
  }

  int i;
  for (i=1; i <= max_trials; ++i)
  {
    system_return_value = system(command.c_str());
    if (system_return_value == 0)
      return i;
    // Else
    // Clean up the unsuccessful run:
    //    sleep(1); // Sleep for 1 second to allow system to close all files.

    //    remove(exonerate_vulgar_filename.c_str());
    //    remove(exonerate_log_filename.c_str());
  }
  return max_trials +1;
}


int run_and_parse_exonerate_for_single_input_file(
          const                          exonerate_parameters &ep,
          vector<vulgar *>               &vec_of_hits_as_in_file,  // New results are appended to this vector
          multimap<faststring, vulgar *> &map_of_vulgar_hits_targets_as_keys, // New results are appended to this map
          vector<stats_for_given_target> &vector_of_hit_stats_for_query_reference) // New results are added to this vec. One entry per reference.
{
  bool     vulgar_file_read_successfully = false;
  unsigned count_appempts_to_read_vulgar_file = 0;

  while (!vulgar_file_read_successfully)
  {
    // Reset vec and map and counters. They might have data from an incomplete vulgar file that we tried to read before.
    //    vec_of_hits_as_in_file.clear();
    //    map_of_vulgar_hits_targets_as_keys.clear();


    faststring input_fasta = ep.fasta_reads;
    faststring basename_of_input_file_without_extension = input_fasta.filename_basename_without_extension();

    faststring vulgar_output_file_name = ep.get_exonerate_vulgar_filename();

    //*****************************
    // Run exonerate if necessary:
    //*****************************

    // Do we need to compute the vulgar file?
    if ( !fileExists(vulgar_output_file_name.c_str()))
    {
      if (global_verbosity >= 1 && count_appempts_to_read_vulgar_file == 1)
        cout << "No vulgar file with the specified name has been found. So the binary \""
        << global_exonerate_binary << "\" will be used to create the vulgar file." << endl;

      int num_trials = run_exonerate(ep);

      if (num_trials == ep.max_trials+1)
      {
        cerr << "ERROR: Running exonerate failed. "
        "The generated vulgar file is incomplete and should be removed manually. Exiting.\n";
        good_bye_and_exit(-1);
      }

      if (global_verbosity >= 1)
        cout << "Exonerate result file was created successfully after N = " << num_trials << " trials." << endl;
    }
    else
    {
      if (global_verbosity >= 1)
        cout << "The specified vulgar file exists, so it does not have to be recomputed with exonerate." << endl;
    }

    //*****************************
    // Parse exonerate output:
    //*****************************

    CFile exonerate_output_file;
    //  map<faststring, vulgar>  m_exonerate_results;

    exonerate_output_file.ffopen(vulgar_output_file_name.c_str());

    if (exonerate_output_file.fail())
    {
      cerr << "ERROR: Could not open file " << vulgar_output_file_name << endl
           << "This file should have been created by exonerate. Its existinace has just been verified. "
           << "Now it cannot be opened any more. Exiting.\n\n";
      good_bye_and_exit(-6);
    }

    faststring line;
    while(true) // Loop over all lines in file and create vulgar objects for all vulgar lines:
    {
      //    char c = exonerate_output_file.getchar();
      //    cout << "char: " << c << endl;

      exonerate_output_file.getline(line);

      if (exonerate_output_file.fail())
        break;

      line.removeSpacesBack();

      // These are valid hits:
      if (line.size() > 7 && line.substr(0,7) == "vulgar:")
      {
        //       cout << "Adding line: " << line << endl;
        // Create vulgar object:
        vulgar *pv = new vulgar(line);

        // Preselection of hits:
        // Depeding on the program parameters, some hits will not be considered. This is a good
        // point to remove them.
        {
          vulgar &vul = *pv;
          short refnumber = pv->get_queryIndex();
          vector_of_hit_stats_for_query_reference[refnumber].increment_hits_in_vulgar_file();

          // Currently, considering frameshift hits is not supported. Anyway we filter them first.
          if (vul.has_F() && !global_include_frameshift_alignments)
          {
            vector_of_hit_stats_for_query_reference[refnumber].increment_skipped_F();
            vector_of_hit_stats_for_query_reference[refnumber].increment_not_considered();
            //      ++skipped_F;
            //      ++count_not_considered;
            if (!vul.has_G() )
              vector_of_hit_stats_for_query_reference[refnumber].increment_skipped_G();
            //        ++skipped_G;
            delete pv;
            continue;
          }

          if ((!vul.has_G() && global_gap_frameshift_mode == 2) )
          {
            vector_of_hit_stats_for_query_reference[refnumber].increment_skipped_no_G();
            //      ++skipped_no_G;
            vector_of_hit_stats_for_query_reference[refnumber].increment_not_considered();
            //      ++count_not_considered;
            delete pv;
            continue;
          }

          if ((vul.has_G() && global_gap_frameshift_mode == 0) )
          {
            vector_of_hit_stats_for_query_reference[refnumber].increment_not_considered();
            //      ++count_not_considered;
            vector_of_hit_stats_for_query_reference[refnumber].increment_skipped_G();
            //      ++skipped_G;
            delete pv;
            continue;
          }

          if (global_relative_score_threshold && vul.relative_score() < global_relative_score_threshold)
          {
            cerr << "NOTE: Exonerate hit skipped due to low relative alignment score: " << vul.relative_score()
            << " " << vul.get_targetID() << endl;
            vector_of_hit_stats_for_query_reference[refnumber].increment_skipped_relative_score();
            vector_of_hit_stats_for_query_reference[refnumber].increment_not_considered();
            //      ++skipped_relative_score;
            //      ++count_not_considered;
            delete pv;
            continue;
          }
        } // END Preselection of hits

        vec_of_hits_as_in_file.push_back(pv); // Removed hits are not added.
        map_of_vulgar_hits_targets_as_keys.emplace(pv->get_targetID(),pv);

        ////        vec_exonerate_results.push_back(make_pair(pv->get_targetID(),pv));
        ////        faststring hit_key = pv->get_hitkey();
        ////        add_or_count(map_exonerate_count_results, hit_key);
      } // END if (line.size() > 7 && line.substr(0,7) == "vulgar:")
      else  // These lines do not contain hits:
      {
        if (line.size() >= 13 && line.substr(0,13) == "Command line:")
          continue;
        if (line.size() >= 9 && line.substr(0,9) == "Hostname:")
          continue;
        if (line.size() >= 31 && line.substr(0,31) ==  "-- completed exonerate analysis")
        {
          vulgar_file_read_successfully = true;
          break;
        }
        if (line.empty()) // Skip empty lines even though they normally do not occur in vulgar files.
          continue;
        cerr << "WARNING: Ignoring non standard line in vulgar file: \"" << line << "\""<< endl;
      }
    } // END  while(true) // Reading all lines of the vulgar file.

    if (!vulgar_file_read_successfully)
    {
      cerr << "WARNING: Vulgar file is incomplete. It does not end with \"-- completed exonerate analysis\"\n";
      cerr << "This is attempt number " << count_appempts_to_read_vulgar_file
      << " to read the vulgar file. I will try to recompute the vulgar file.\n";
      exonerate_output_file.ffclose();
      remove(vulgar_output_file_name.c_str());
    }
    else
    {
      exonerate_output_file.ffclose();
    }

    ++count_appempts_to_read_vulgar_file;
  } // END  while (!vulgar_file_read_successfully)

  // Debug output:
  if (/* DISABLES CODE */ (0))
  {
    cout << "Finished reading exonerate results. Please hit any key to continue.\n" << endl;
    cin.get();
  }

  return 0;
}
