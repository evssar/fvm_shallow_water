/*

Function for reading .grid file

*/

#include "sw.h"

void read_grid(sw &sw)
{

    std::string file_line;
    std::smatch file_line_match;
    std::regex file_line_regex;

    sw.file.open(sw.case_name + ".grid", std::fstream::in);

    if (sw.file.is_open())
    {
        file_line_regex.assign("^([0-9]+),([0-9]+)");
        std::getline(sw.file, file_line);

        if (std::regex_search(file_line, file_line_match, file_line_regex))
        {
            int N_vertices, N_faces;
            double x, y, z;
            Eigen::MatrixXd vertices;
            Eigen::MatrixXi faces;
            N_vertices = std::stoi(file_line_match[1]);
            N_faces = std::stoi(file_line_match[2]);

            vertices = Eigen::MatrixXd::Zero(N_vertices, 3);
            faces = Eigen::MatrixXi::Zero(N_faces, 4);

            

            file_line_regex.assign("^(-?[0-9]+\\.[0-9]+),(-?[0-9]+\\.[0-9]+),(-?[0-9]+\\.[0-9]+)");

            for (unsigned i = 0; i < N_vertices; i++)
            {
                std::getline(sw.file, file_line);

                if (!file_line.empty())
                {
                    if (std::regex_search(file_line, file_line_match, file_line_regex))
                    {
                        vertices(i, 0) = std::stod(file_line_match[1]);
                        vertices(i, 1) = std::stod(file_line_match[2]);
                        vertices(i, 2) = std::stod(file_line_match[3]);
                    }
                }
            }

    
            file_line_regex.assign("^([0-9]+),([0-9]+),([0-9]+),([0-9]+)");
            

            for (unsigned i = 0; i < N_faces; i++)
            {
                std::getline(sw.file, file_line);

                if (!file_line.empty())
                {
                    if (std::regex_search(file_line, file_line_match, file_line_regex))
                    {
                        faces(i, 0) = std::stoi(file_line_match[1]);
                        faces(i, 1) = std::stoi(file_line_match[2]);
                        faces(i, 2) = std::stoi(file_line_match[3]);
                        faces(i, 3) = std::stoi(file_line_match[4]);
                    }
                }
            }

            Eigen::Vector3d r1, r2, r3, r4;
            for (int i = 0; i < N_faces; i++)
            {

                    r1(0) = vertices(faces(i, 0), 0);
                    r1(1) = vertices(faces(i, 0), 1);
                    r1(2) = vertices(faces(i, 0), 2);

                    r2(0) = vertices(faces(i, 1), 0);
                    r2(1) = vertices(faces(i, 1), 1);
                    r2(2) = vertices(faces(i, 1), 2);

                    r3(0) = vertices(faces(i, 2), 0);
                    r3(1) = vertices(faces(i, 2), 1);
                    r3(2) = vertices(faces(i, 2), 2);

                    r4(0) = vertices(faces(i, 3), 0);
                    r4(1) = vertices(faces(i, 3), 1);
                    r4(2) = vertices(faces(i, 3), 2);

                    cell cell(r1, r2, r3, r4);

                   
                    sw.cells.push_back(cell);
                
            }
        }

        sw.file.close();

        sw.N_cells = sw.cells.size();
    }
    else
    {
        std::cout << "Failed reading grid file!" << std::endl;
        std::exit(-1);
    }
}
