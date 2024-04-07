#include "hierarchical_overlap_graph.hpp"

/*
TODO LIST:
 - passare da std::string a struttura più efficiente con 2 bit per carattere
 - usare bucketsort per ordinare le stringhe in ordine di lunghezza decresente
 - controllare se è possibile non fare uso delle strutture R_l
*/

/* Constructor of class HOG */
HOG::HOG(vector<string> &in_reads)
{
    // Import the reads inside HOG
    reads.swap(in_reads);

    // Creation of root node
    parent.push_back(-1);
    sibling.push_back(-1);
    read_id.push_back(0);
    read_length.push_back(0);
    leaves.push_back(true);
    _size = 1;

    auto start = std::chrono::high_resolution_clock::now();
    std::sort(reads.begin(), reads.end());
    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "Time to sort reads: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << "ms" << std::endl;

    auto comp_start = std::chrono::high_resolution_clock::now();

    start = std::chrono::high_resolution_clock::now();
    this->add_reads();
    std::cout << "Number of nodes in the Trie: " << _size << std::endl;
    stop = std::chrono::high_resolution_clock::now();
    std::cout << "Time to insert reads: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << "ms" << std::endl;
    memory_usage_MB("Memory occupation after Trie creation: ");

    start = std::chrono::high_resolution_clock::now();
    this->add_failure_links();
    stop = std::chrono::high_resolution_clock::now();
    std::cout << "Time to add failure_links: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << "ms" << std::endl;
    memory_usage_MB("Memory occupation after failure_link creation: ");

    start = std::chrono::high_resolution_clock::now();
    this->mark_ehog();
    stop = std::chrono::high_resolution_clock::now();
    std::cout << "Time to mark EHOG nodes: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << "ms" << std::endl;
    memory_usage_MB("Memory occupation after EHOG marking: ");

    start = std::chrono::high_resolution_clock::now();
    this->contract();
    stop = std::chrono::high_resolution_clock::now();
    std::cout << "Time to contract EHOG: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << "ms" << std::endl;
    memory_usage_MB("Memory occupation after EHOG contraption: ");

    start = std::chrono::high_resolution_clock::now();
    this->mark_Rl();
    stop = std::chrono::high_resolution_clock::now();
    std::cout << "Time to create R_l lists: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << "ms" << std::endl;
    memory_usage_MB("Memory occupation after creation of R_l lists: ");

    start = std::chrono::high_resolution_clock::now();
    this->mark_hog();
    stop = std::chrono::high_resolution_clock::now();
    std::cout << "Time to mark EHOG nodes: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << "ms" << std::endl;
    memory_usage_MB("Memory occupation after HOG marking: ");

    start = std::chrono::high_resolution_clock::now();
    this->contract();
    stop = std::chrono::high_resolution_clock::now();
    std::cout << "Time to contract HOG: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << "ms" << std::endl;
    memory_usage_MB("Memory occupation after HOG contraption: ");

    auto comp_stop = std::chrono::high_resolution_clock::now();
    std::cout << "\n\nTotal HOG creation time: " << std::chrono::duration_cast<std::chrono::milliseconds>(comp_stop - comp_start).count() << "ms" << std::endl;
}

/* Returns the index of the first child of the node passed as argument */
uint32_t HOG::get_first_child(uint32_t node)
{
    if (leaves[node])
    {
        return -1;
    }
    return node + 1;
}

/* Returns the index of the last child of the node passed as argument */
uint32_t HOG::get_last_child(uint32_t node)
{
    if (leaves[node])
    {
        return -1;
    }

    uint32_t child = node + 1;
    while (sibling[child] != -1)
    {
        child = sibling[child];
    }

    return child;
}

/* Returns the index of the child of the node passed as argument that represents the char c or -1 if that child doesn't exist */
uint32_t HOG::get_child(uint32_t node, char c)
{
    if (leaves[node])
    {
        return -1;
    }

    uint32_t child = node + 1;
    while (child != -1 && reads[read_id[child]][read_length[child] - 1] != c)
    {
        child = sibling[child];
    }

    return child;
}

/* Returns the string label that the node represents */
string HOG::get_label(uint32_t node)
{
    string label = "";

    if (node != 0)
    {
        uint32_t from = read_length[parent[node]];
        uint32_t length = read_length[node] - read_length[parent[node]];

        label = reads[read_id[node]].substr(from, length);
    }

    return label;
}

/* Adds the reads in the Trie-like structure */
void HOG::add_reads()
{
    for (uint32_t i = 0; i < reads.size(); i++)
    {
        string &read = reads[i];
        uint32_t p = 0; // parent index starts from root every time

        for (size_t j = 0; j < read.size(); j++)
        {
            uint32_t c = get_child(p, read[j]);

            if (c == -1) // doesn't have child of label ch
            {
                c = _size;

                parent.push_back(p);
                sibling.push_back(-1);

                read_id.push_back(i);
                read_length.push_back(j + 1);

                leaves.push_back(true);

                leaves[p] = false;

                uint32_t prev_sibling = get_last_child(p);
                if (prev_sibling != c) // current is not only first child
                {
                    sibling[prev_sibling] = c;
                }

                _size++;
            }

            p = c;
        }
    }
}

/* Adds the failure links to the nodes */
void HOG::add_failure_links()
{
    failure_link = vector<uint32_t>(_size);

    queue<uint32_t> q;
    q.push(0);

    while (q.size())
    {
        uint32_t curr = q.front(); // Estraggo indice del primo nodo dalla coda
        q.pop();

        while (curr != -1)
        {
            // Aggiungo alla coda il figlo del nodo se c'e'
            uint32_t c = get_first_child(curr);
            if (c != -1)
            {
                q.push(c);
            }

            // Se nodo radice, fail_link punta -1 (base case)
            if (parent[curr] == -1)
            {
                failure_link[curr] = -1;
                curr = sibling[curr];
                continue;
            }
            // Se nodo figlio della radice, fail_link punta alla radice (base case)
            if (parent[curr] == 0)
            {
                failure_link[curr] = 0;
                curr = sibling[curr];
                continue;
            }

            uint32_t parent_fail = failure_link[parent[curr]]; // Index of next longest suffix node of parent

            // Risalgo la catena di failure link dal padre finche' non trovo
            // un nodo che ha come figlio la stessa label o non trovo la radice

            char current_label = reads[read_id[curr]][read_length[curr] - 1]; // Non e' un problema il fatto che la radice non ha caretteri perche' la radice e' un caso base
            while (get_child(parent_fail, current_label) == -1 && parent_fail != 0)
            {
                parent_fail = failure_link[parent_fail];
            }

            if (parent_fail == 0 && get_child(0, current_label) == -1)
            {
                failure_link[curr] = 0;
            }
            else
            {
                failure_link[curr] = get_child(parent_fail, current_label);
            }

            curr = sibling[curr];
        }
    }
}

/* Marks the only nodes that will remain to form the EHOG */
void HOG::mark_ehog()
{
    marked = vector<bool>(_size, false);

    for (uint32_t i = 0; i < _size; i++)
    {
        if (leaves[i])
        {
            uint32_t ref = i;

            while (ref != -1)
            {
                if (marked[ref])
                {
                    break;
                }

                marked[ref] = true;

                ref = failure_link[ref];
            }
        }
    }
}

/* Contracts the Trie to eliminate the unmarked nodes and merge the nodes that need to merge with their ancestors */
void HOG::contract()
{
    vector<uint32_t> new_parent{};
    vector<uint32_t> new_sibling{};
    vector<uint32_t> new_failure_link{};
    vector<bool> new_leaves{};

    vector<uint32_t> new_read_id{};
    vector<uint32_t> new_read_length{};

    uint32_t new_size = 0;

    vector<uint32_t> failure_link_conversion(_size);

    std::function<void(uint32_t, uint32_t)> contract_subtree = [&](uint32_t local_root, uint32_t parent_of_local_root)
    {
        stack<uint32_t> local_to_visit;
        local_to_visit.push(local_root);

        while (local_to_visit.size())
        {
            uint32_t n = local_to_visit.top();
            local_to_visit.pop();

            while (!marked[n])
            {
                n = get_first_child(n); // no need to check if is -1 (parent is a leave), because leaves are always marked

                uint32_t s = sibling[n];
                while (s != -1)
                {
                    local_to_visit.push(s);
                    s = sibling[s];
                }
            }

            new_parent.push_back(parent_of_local_root);
            new_sibling.push_back(-1);

            uint32_t last_child = parent_of_local_root + 1;
            while (new_sibling[last_child] != -1)
            {
                last_child = new_sibling[last_child];
            }
            if (last_child != new_size) // current is not only first child
            {
                new_sibling[last_child] = new_size;
            }

            new_leaves.push_back(leaves[n]);
            new_read_id.push_back(read_id[n]);
            new_read_length.push_back(read_length[n]);

            uint32_t f = failure_link[n];
            while (f != -1 && !marked[f])
            {
                f = failure_link[f];
            }
            new_failure_link.push_back(f);
            failure_link_conversion[n] = new_size;

            uint32_t current_size = new_size;
            new_size++;

            uint32_t ch = get_first_child(n);
            while (ch != -1)
            {
                contract_subtree(ch, current_size);
                ch = sibling[ch];
            }
        }
    };

    contract_subtree(0, -1);

    for (uint32_t i = 0; i < new_failure_link.size(); i++)
    {
        if (new_failure_link[i] != -1)
        {
            new_failure_link[i] = failure_link_conversion[new_failure_link[i]];
        }
    }

    parent.swap(new_parent);
    sibling.swap(new_sibling);
    failure_link.swap(new_failure_link);
    leaves.swap(new_leaves);

    read_id.swap(new_read_id);
    read_length.swap(new_read_length);

    _size = new_size;

    std::cout << "Trie node number after contraption: " << _size << std::endl;
}

/* Creates the R_l lists for all nodes */
void HOG::mark_Rl()
{
    R_l = vector<vector<uint32_t>>(_size, vector<uint32_t>());

    for (uint32_t i = 0; i < _size; i++)
    {
        if (leaves[i])
        {
            uint32_t failure = failure_link[i];

            while (failure != -1)
            {
                R_l[failure].push_back(i);
                failure = failure_link[failure];
            }
        }
    }
}

/* Marks the only nodes that will remain to form the HOG */
void HOG::mark_hog()
{
    marked = vector<bool>(_size, false);
    for (uint32_t i = 0; i < _size; i++)
    {
        marked[i] = leaves[i];
    }
    marked[0] = true;

    vector<bool> inS = vector<bool>(_size, false);
    set<uint32_t> S;
    unordered_map<uint32_t, stack<uint32_t>> Sx;
    for (uint32_t i = 0; i < _size; i++)
    {
        if (leaves[i])
        {
            Sx[i] = stack<uint32_t>();
        }
    }

    vector<bool> visited(_size, false); // to distinguish between first visit and second (last) visit
    stack<uint32_t> to_visit;
    to_visit.push(0);

    while (to_visit.size())
    {
        uint32_t v = to_visit.top();

        if (visited[v] || leaves[v])
        {
            to_visit.pop();
        }

        if (leaves[v]) // v is a leave
        {
            for (uint32_t x : S)
            {
                marked[Sx[x].top()] = true;
                inS[x] = false;
            }

            S.clear();
        }
        else // v is an internal node
        {
            if (!visited[v]) // First visit
            {
                for (uint32_t x : R_l[v])
                {
                    Sx[x].push(v);

                    if (inS[x] == false)
                    {
                        inS[x] = true;
                        S.insert(x);
                    }
                }

                // Visit in DFS
                uint32_t c = get_first_child(v);
                while (c != -1)
                {
                    to_visit.push(c);
                    c = sibling[c];
                }
            }
            else // Last visit
            {
                for (uint32_t x : R_l[v])
                {
                    Sx[x].pop();
                    if (Sx[x].size() != 0 && marked[Sx[x].top()] == false)
                    {
                        if (inS[x] == false)
                        {
                            inS[x] = true;
                            S.insert(x);
                        }
                    }
                    else
                    {
                        if (inS[x] == true)
                        {
                            inS[x] = false;
                            S.erase(x);
                        }
                    }
                }
            }
        }

        visited[v] = true;
    }

    vector<vector<uint32_t>>().swap(R_l); // better free R_l memory since not useful after this, and next contraction will require additional memory
}

/* Prints in standars output the nodes in order of index */
void HOG::print()
{
    std::cout << std::endl
              << "size: " << _size
              << std::endl;

    queue<uint32_t> q;
    q.push(0);

    while (q.size())
    {
        uint32_t curr = q.front();
        q.pop();

        while (curr != -1)
        {
            q.push(get_first_child(curr));
            print_node(curr);
            curr = sibling[curr];
        }
    }
}

/* Prints the node's current status */
void HOG::print_node(uint32_t i)
{

    std::cout << std::endl
              << "node_number:" << i << ' '
              << "read_id:" << read_id[i] << ' '
              << "read_len:" << read_length[i] << ' '
              << "string_label:" << get_label(i) << ' '
              << "parent:" << parent[i] << ' '
              << "sibling:" << sibling[i] << ' '
              << "failure_link:" << failure_link[i] << ' '
              << "marked:" << marked[i] << ' '
              << "leave:" << leaves[i] << ' ';
    // cout << "R_l:";
    // for (uint32_t x : R_l[i])
    //     cout << x << ' ';

    std::cout << std::endl;
}
