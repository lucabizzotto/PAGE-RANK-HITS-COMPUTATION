#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <map>
#include <sys/mman.h>
#include <cstdlib>
#include <cmath>


typedef std::pair<unsigned int, unsigned int> pair;


class InDegree
{
public:
    InDegree(std::string name_file) : NAME_FILE(name_file)
    {
        // determine the number of nodes and edges in the graph
        number_of_nodes_and_edges();
        // initialize vector
        this->initialize_in_Deg();
        // read file and save values
        this->initialize();
        std::stable_sort(this->ptr, this->ptr + this->N_edges, Compararator_first_element());
        std::stable_sort(this->ptr, this->ptr + this->N_edges, Compararator_second_element());
    }


    // vector of prestige
    std::vector<unsigned int> in_Deg;
    // number of interations to converge
    unsigned int n_convergence = 0;


    // custom class to execute a custom sort
    class Compararator_first_element {
    public:
        bool operator()(const pair& a, const pair& b) {
            return a.first < b.first;
        }
    };


    // custom class to execute a custom sort
    class Compararator_second_element {
    public:
        bool operator()(const pair& a, const pair& b) {
            return a.second < b.second;
        }
    };


    // compute pagerank
    void compute()
    {
        for (unsigned int i = 0; i < this->N_edges; i++)
        {
            in_Deg[this->ptr[i].second - this->Min_node] += 1;
        }
    }


private:

    // holds the path of file
    const std::string NAME_FILE;
    // pointer to permanent memory
    pair* ptr;
    // holds the dimension of ptr and the number of edges
    unsigned int N_edges;
    // holds the number of nodes
    unsigned int N_nodes;
    // holds the value of min node
    unsigned int Min_node;
    // holds the value of max node
    unsigned int Max_node;


    // initialize p_k at time 0
    inline void initialize_in_Deg()
    {
        this->in_Deg.clear();
        this->in_Deg.resize(this->N_nodes, 0);
    }


    // open file
    inline std::ifstream open_file()
    {
        // create instance of file
        std::ifstream file;
        // open file
        file.open(this->NAME_FILE);
        // control the validity of file
        if (!file.is_open())
        {
            throw std::runtime_error("Could not open file");
        }
        return file;
    }


    // close file
    inline void close_file(std::ifstream& file)
    {
        file.close();
    }


    // return the number of both nodes and edges in the graph
    void number_of_nodes_and_edges()
    {
        std::ifstream file = this->open_file();
        // create string to read file
        std::string line;
        // create string to store all info of file
        std::string info_string;
        // read file
        while (std::getline(file, line))
        {   // in this way we read only the description of the file
            if (line[0] != '#')
            {
                break;
            }
            else
            {
                // append all lines of description in a single string
                info_string.append(line);
            }
        }
        // close the file
        this->close_file(file);
        // create rule for regexp to find edges
        std::regex r_edges("(Edges: [0-9]*)");
        // create rule for regexp to find nodes
        std::regex r_nodes("(Nodes: [0-9]*)");
        std::regex r_min_node("(Min_node: [0-9]*)");
        std::regex r_max_node("(Max_node: [0-9]*)");
        // find the numbers
        this->N_edges = regex_exec(r_edges, info_string);
        this->N_nodes = regex_exec(r_nodes, info_string);
        this->Min_node = regex_exec(r_min_node, info_string);
        this->Max_node = regex_exec(r_max_node, info_string);
    }


    // return the value associate with regex and string
    inline int regex_exec(std::regex r, std::string info_string)
    {
        // array for regexp result
        std::smatch m;
        // compute the regexp on text
        std::regex_search(info_string, m, r);
        // retrieve the result of regexp
        std::string result = m.str(0);
        // split the result in strings with a specific break in this way we can take only the number and transform it in an integer
        std::istringstream result_line(result);
        std::vector<std::string> v_result;
        for (std::string piece; std::getline(result_line, piece, ' ');)
        {
            v_result.push_back(piece);
        }
        // return the number
        return std::stoi(v_result[1]);
    }


    // from the line return the pair that compose the edge
    inline pair take_edge_from_line(std::string& line)
    {
        // split the result in strings with a specific break
        std::istringstream result_line(line);
        // contains the values for the pair
        std::vector<int> v_result;
        for (std::string piece; std::getline(result_line, piece, '\t');)
        {
            v_result.push_back(std::stoi(piece));
        }
        return pair(v_result[0], v_result[1]);
    }


    // read the file and save in permanent memory
    void initialize()
    {
        // allocate the right amount of space in permanent memory to save the edges of graph as a set of pairs
        this->ptr = (pair*)mmap(NULL, this->N_edges * sizeof(pair), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0, 0);
        // control if the allocate was done in right way
        if (ptr == MAP_FAILED)
        {
            throw std::runtime_error("Mapping Failed\n");
        }
        // create instance of file
        std::ifstream file = this->open_file();
        // create string to read each row of file
        std::string line;
        // take the position of the ptr to insert the pair
        unsigned int i = 0;
        // read file
        while (std::getline(file, line))
        {   // in this way we don't read the description of the file
            if (line[0] != '#')
            {
                // store the pair that define the edge
                this->ptr[i] = this->take_edge_from_line(line);
                // increment the pointer
                i++;
            }
        }
        // close file
        this->close_file(file);
    }

};