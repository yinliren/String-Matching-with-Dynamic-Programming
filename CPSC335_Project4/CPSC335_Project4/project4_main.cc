
#include "project4.hh"
#include "timer.hh"

int main() {
    ProteinVector proteins;
    auto load_successful = load_proteins(proteins, "proteins_large.txt");
    assert( load_successful );

    BlosumPenaltyArray bpa;
    load_successful = load_blosum_file(bpa, "blosum62.txt");
    assert( load_successful );

    std::vector<std::string> testProteins;
    testProteins.push_back("PIEPCMGA");
    testProteins.push_back("TQGASNIGE");
    testProteins.push_back("ALAKLIRYGG");
    testProteins.push_back("CSNPNLSDFGR");
    testProteins.push_back("MYPEPTIDE");

    std::cout << "------------------- Dynamic Programming ------------------" << std::endl;
    for (int i = 0; i < testProteins.size(); i++) {
        std::string align_string1;
        std::string align_string2;
        int best_Score;
        std::string searchString =     testProteins[i];
        Timer timer;
        std::cout << "String to Match = " << testProteins[i] << std::endl;
        std::shared_ptr<Protein> best_protein = local_alignment_best_match(proteins,
                                                                          searchString,
                                                                          bpa,
                                                                          align_string1,
                                                                          align_string2,
                                                                           best_Score);
        std::cout << "Best Protein Description: " << best_protein->description << std::endl;
        std::cout << align_string1 << std::endl;
        std::cout << align_string2 << std::endl;
        std::cout << "Best Score: " << best_Score << std::endl;
        std::cout << "Time Elapsed: " << timer.elapsed() << std::endl;
    }

//    //=====================================UNIT TEST===================================================
//
//    ProteinVector proteins;
//    auto load_successful = load_proteins(proteins, "proteins_large.txt");
//    assert( load_successful );
//
//    BlosumPenaltyArray bpa;
//    load_successful = load_blosum_file(bpa, "blosum62.txt");
//    assert( load_successful );
//
////    std::vector<std::string> testProteins;
////    testProteins.push_back("KSSTTSSTSKKPA");
//
//    std::string align_string1;
//    std::string align_string2;
//    std::string searchString = "RLIKKTTGSSSSSSSKKKDKEKEKEKSSTTSKKPASASSSSHGTTHSSASSTGSKSTTEKGKQSGSVPSQGKHHSSSTSKTKTATTPSSSSSSSRSSSVSRSGSSSTKKTSSRKGQEQSKQSQQPSQSQKQGSSSSSAAIMNPTP";
//    Timer timer;
//    std::cout << "String to Match = " << searchString << std::endl;
////    int best_score = local_alignment("MMRGFKQRLIKK",
////                    searchString,
////                    bpa,
////                    align_string1,
////                    align_string2);
//    std::shared_ptr<Protein> best_protein = local_alignment_best_match(proteins, searchString,
//                               bpa,
//                               align_string1,
//                               align_string2);
//
//    std::cout << best_protein->description << std::endl;
//    std::cout << align_string1 << std::endl;
//    std::cout << align_string2 << std::endl;
////    std::cout << best_score << std::endl;
//    std::cout << timer.elapsed() << std::endl;
//
//

//    >sp|P38903|2A5D_YEAST Serine/threonine-protein phosphatase 2A 56 kDa regulatory subunit delta isoform OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) GN=RTS1 PE=1 SV=2 MMRGFKQRLIKKTTGSSSSSSSKKKDKEKEKEKSSTTSSTSKKPASASSSSHGTTHSSASSTGSKSTTEKGKQSGSVPSQGKHHSSSTSKTKTATTPSSSSSSSRSSSVSRSGSSSTKKTSSRKGQEQSKQSQQPSQSQKQGSSSSSAAIMNPTPVLTVTKDDKSTSGEDHAHPTLLGAVSAVPSSPISNASGTAVSSDVENGNSNNNNMNINTSNTQDANHASSQSIDIPRSSHSFERLPTPTKLNPDTDLELIKTPQRHSSSRFEPSRYTPLTKLPNFNEVSPEERIPLFIAKVDQCNTMFDFNDPSFDIQGKEIKRSTLDELIEFLVTNRFTYTNEMYAHVVNMFKINLFRPIPPPVNPVGDIYDPDEDEPVNELAWPHMQAVYEFFLRFVESPDFNHQIAKQYIDQDFILKLLELFDSEDIRERDCLKTTLHRIYGKFLSLRSFIRRSMNNIFLQFIYETEKFNGVAELLEILGSIINGFALPLKEEHKVFLVRILIPLHKVRCLSLYHPQLAYCIVQFLEKDPLLTEEVVMGLLRYWPKINSTKEIMFLNEIEDIFEVIEPLEFIKVEVPLFVQLAKCISSPHFQVAEKVLSYWNNEYFLNLCIENAEVILPIIFPALYELTSQLELDTANGEDSISDPYMLVEQAINSGSWNRAIHAMAFKALKIFLETNPVLYENCNALYLSSVKETQQRKVQREENWSKLEEYVKNLRINNDKDQYTIKNPELRNSFNTASENNTLNEENENDCDSEIQ
    
    
   // >sp|P48016|ATH1_YEAST Vacuolar acid trehalase OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) GN=ATH1 PE=1 SV=1
    //MKRIRSLWFNAEASYSNLNNSPSLRNKNSTGNNSRSKNYRSFSRFDLINSILLLMMLFLLAIFVTALYLTKSSRLTYSHASRAALFNPLGVISPSLGNHTLNYDPEARESSKKLYELLSDFNTAYYDDENMILGSNLFSKNTYSRQPYVANGYIGSRIPNIGFGYALDTLNFYTDAPGALNNGWPLRNHRFAGAFVSDFYCLQPKLNSTNFPELDDVGYSTVISSIPQWTNLQFSLVNDSKWFNPQNVTLDDVTNYSQNLSMKDGIVTTELDWLNSQIHVKSEIWAHRHIHPLGVVSLEISLNTDHLPSDFDSLDVNIWDILDFNTSHRTVLHSTGTDEKNNAVFMIVQPDNVPSSNCAIYSTCTVKYENSTNPINSSESFEEKDVSSNIYNVILTEDQPKIIVHKYVGIMSTEFNKNKEQQDNTNIGLAKMIALNSKGNYEKLLSSHKRAWYDLYNDAFIEIPSDSLLEMTARSSLFHLLANTRDYNVSSDRGLPVGVSGLSSDSYGGMVFWDADIWMEPALLPFFPNVAQNMNNYRNATHSQAKLNAEKYGYPGAIYPWTSGKYANCTSTGPCVDYEYHINVDVAMASFSIYLNGHEGIDDEYLRYTTWPIIKNAAQFFTAYVKYNSSLGLYETYNLTDPDEFANHINNGAFTNAGIKTLLKWATDIGNHLGEVVDPKWSEISKDIYIPRSSSNITLEYSGMNSSVEIKQADVTLMVYPLGYINDESILNNAIKDLYYYSERQSASGPAMTYPVFVAAAAGLLNHGSSSQSYLYKSVLPYLRAPFAQFSEQSDDNFLTNGLTQPAFPFLTANGGFLQSILFGLTGIRYSYEVDPDTKKINRLLRFNPIELPLLPGGIAIRNFKYMNQVLDIIIDDHNGTIVHKSGDVPIHIKIPNRSLIHDQDINFYNGSENERKPNLERRDVDRVGDPMRMDRYGTYYLLKPKQELTVQLFKPGLNARNNIAENKQITNLTAGVPGDVAFSALDGNNYTHWQPLDKIHRAKLLIDLGEYNEKEITKGMILWGQRPAKNISISILPHSEKVENLFANVTEIMQNSGNDQLLNETIGQLLDNAGIPVENVIDFDGIEQEDDESLDDVQALLHWKKEDLAKLIEQIPRLNFLKRKFVKILDNVPVSPSEPYYEASRNQSLIEILPSNRTTFTIDYDKLQVGDKGNTDWRKTRYIVVAVQGVYDDYDDDNKGATIKEIVLND

  return 0;
}




