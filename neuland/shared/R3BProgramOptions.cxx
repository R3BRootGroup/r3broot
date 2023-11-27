#include "R3BProgramOptions.h"

namespace R3B
{
    auto ProgramOptions::Verify(int argc, const char** argv) -> bool
    {

        try
        {
            po::store(po::command_line_parser(argc, argv).positional(pos_desc_).options(desc_).run(), varMap_);
            po::notify(varMap_);
            for (auto& registrie : registries_)
            {
                registrie.second->Retrieve(varMap_);
            }
        }
        catch (const std::runtime_error& err)
        {
            std::cerr << "An exception occurs: " << err.what() << std::endl;
            return false;
        }
        catch (const std::exception& err)
        {
            // std::cerr << "exception occurs: " << boost::diagnostic_information(err) << std::endl;
            std::cerr << "An exception occurs: " << err.what() << std::endl;
            return false;
        }
        catch (...)
        {
            std::cerr << "An unrecognizable exception occurs!" << std::endl;
            return false;
        }
        return true;
    }
} // namespace R3B
