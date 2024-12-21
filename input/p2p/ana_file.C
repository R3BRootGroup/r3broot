
void ana_file()
{

    std::ifstream fFile; //! Input file handle

    const std::string fFileName = "quasi_ascii.txt"; //! Input file name

    // Event variable to be read from file
    int eventId = -1;
    int nTracks = -1;
    // Track variables to be read from file
    int iPid = -1;
    int iZ = 0;
    int iA = 0;
    double px = 0.;
    double py = 0.;
    double pz = 0.;
    double vx = 0.;
    double vy = 0.;
    double vz = 0.;
    double ivx = 0.;
    double ivy = 0.;
    double ivz = 0.;

    auto fheang_1 = new TH2F("fheang_1", "", 200, 0, 90, 2000, 0, 800);

    if (!fFile.is_open())
    {

        fFile.open(fFileName);
    }

    while (!fFile.eof())
    {
        // Read event number and number of primary particles.
        if (!(fFile >> eventId >> nTracks))
        {
        }
        // Ignore the other stuff that might still be on that line
        fFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        // Loop over tracks in the current event
        for (int itrack = 0; itrack < nTracks; itrack++)
        {
            if (!(fFile >> iPid >> iZ >> iA >> px >> py >> pz >> ivx >> ivy >> ivz))
            {
            }
            // Ignore the other stuff that might still be on that line
            fFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            // std::cout<< pz<<std::endl;

            fheang_1->Fill(180. * atan(sqrt(px * px + py * py) / pz) / 3.1415,
                           1000. * (sqrt(px * px + py * py + pz * pz + 0.93827 * 0.93827) - 0.93827));

        } //! tracks

        printf("Processed: \033[32m %d \033[0m \033[31m \033[0m \r", eventId);
        fflush(stdout);
    }

    TCanvas* c4 = new TCanvas("c4");
    fheang_1->Draw("col");
}
