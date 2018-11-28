#include "Player.hh"
#include <limits>
#include <queue>
#include <set>
#include <utility>

#define PLAYER_NAME Linearithmic

using namespace std;

class PLAYER_NAME : public Player {
private:
    set<Pos> water_cells;
    set<Pos> station_cells;

    const int infinity = numeric_limits<int>::max();

    bool calculated_nearest_water;
    vector<vector<int>> _nearest_water;
    bool calculated_nearest_station;
    vector<vector<int>> _nearest_station;

    /** A vector containing the locked cells.
     * A lock is valid until the round marked by this vector (not included).
     * A lock may not be breached, since doing so may result in the death of
     * either of the units (or both!).
     */
    vector<vector<int>> cell_locks;

    inline bool lock(int x, int y) {
        if (cell_locks[x][y] <= round()) {
            cell_locks[x][y] = round() + 1;
            return true;
        } else {
            return false;
        }
    }

    inline bool lock(Pos p) { return lock(p.i, p.j); }

    inline bool is_locked(int x, int y) { return cell_locks[x][y] > round(); }

    inline bool is_locked(Pos p) { return is_locked(p.i, p.j); }

    // TODO: Template this
    void calculate_nearest_water() {
        queue<Pos> bfsq;
        for (auto pos : water_cells) {
            _nearest_water[pos.i][pos.j] = 0;
            bfsq.push(pos);
        }
        while (not bfsq.empty()) {
            Pos curr = bfsq.front();
            bfsq.pop();
            int new_dist = _nearest_water[curr.i][curr.j] + 1;
            for (int d = 0; d < DirSize; ++d) {
                Pos next_pos = curr + Dir(d);
                if (not pos_ok(next_pos)) { continue; }
                switch (cell(next_pos).type) {
                    case Station:
                    case Wall:
                        // Distance = inf (impassable)
                        continue;
                    default:
                        break;
                }
                if (_nearest_water[next_pos.i][next_pos.j] > new_dist) {
                    _nearest_water[next_pos.i][next_pos.j] = new_dist;
                    bfsq.push(next_pos);
                }
            }
        }
    }

    void calculate_nearest_station() {
        queue<Pos> bfsq;
        for (auto pos : station_cells) {
            _nearest_station[pos.i][pos.j] = 0;
            bfsq.push(pos);
        }
        while (not bfsq.empty()) {
            Pos curr = bfsq.front();
            bfsq.pop();
            int base_dist = _nearest_station[curr.i][curr.j] + 1;
            for (int d = 0; d < DirSize - 1; ++d) { // Ignore None
                int new_dist = base_dist + 1;
                Pos next_pos = curr + Dir(d);
                if (not pos_ok(next_pos)) { continue; }
                switch (cell(next_pos).type) {
                    case City:
                    case Wall:
                    case Water:
                        // Distance = inf (impassable)
                        continue;
                    case Desert:
                        new_dist += 4; // Penalty for not going through the road
                        break;
                    default:
                        break;
                }
                if (_nearest_station[next_pos.i][next_pos.j] > new_dist) {
                    _nearest_station[next_pos.i][next_pos.j] = new_dist;
                    bfsq.push(next_pos);
                }
            }
        }
    }

    Dir nearest_water(Pos from) {
        if (not calculated_nearest_water) {
            calculate_nearest_water();
            calculated_nearest_water = true;
        }
        int current_min = _nearest_water[from.i][from.j];
        Dir current_dir = None;
        for (int d = 0; d < DirSize - 1; ++d) { // Ignore None
            Pos next_pos = from + Dir(d);
            if (not pos_ok(next_pos)) { continue; }
            int next_pos_dist = _nearest_water[next_pos.i][next_pos.j];
            if (next_pos_dist < current_min) {
                current_min = next_pos_dist;
                current_dir = Dir(d);
            }
        }
        return current_dir;
    }

    static void print_matrix(const vector<vector<int>> &m) {
        for (auto v : m) {
            for (auto i : v) { cerr << i << ' '; }
            cerr << endl;
        }
    }

    Dir nearest_station(Pos from) {
        if (not calculated_nearest_station) {
            calculate_nearest_station();
            calculated_nearest_station = true;
        }
        int current_min = _nearest_station[from.i][from.j];
        Dir current_dir = None;
        vector<int> distances(DirSize - 1);
        for (int d = 0; d < DirSize - 1; ++d) { // Ignore None
            Pos next_pos = from + Dir(d);
            if (not pos_ok(next_pos)) {
                distances[d] = -1;
                continue;
            }
            int next_pos_dist = _nearest_station[next_pos.i][next_pos.j];
            distances[d] = next_pos_dist;
            if (next_pos_dist < current_min) {
                current_min = next_pos_dist;
                current_dir = Dir(d);
            }
        }
        if (false) {
            CellType from_t = cell(from).type;
            Pos to = from + current_dir;
            CellType to_t = cell(to).type;
            if (from_t == Road and to_t == Desert) {
                cerr << "Going from " << from << " [Road] ";
                cerr << "to " << to << " [Desert]" << endl;
                cerr << "Alternatives were:" << endl;
                for (int d = 0; d < DirSize - 1; ++d) {
                    Pos next_pos = from + Dir(d);
                    CellType type = cell(next_pos).type;
                    string type_s;
                    switch (type) {
                        case Road:
                            type_s = "Road";
                            break;
                        case Desert:
                            type_s = "Desert";
                            break;
                        default:
                            type_s = "Invalid";
                            break;
                    }
                    cerr << "  " << next_pos << " [" << type_s
                         << "]: " << distances[d] << endl;
                }
                cin.get();
            }
        }
        return current_dir;
    }

    void move_warrior(int id) {
        Pos pos = unit(id).pos;
        Dir dir = nearest_water(pos);
        pos += dir;
        if (lock(pos)) { command(id, dir); }
    }

    void move_car(int id) {
        Pos pos = unit(id).pos;
        Dir dir = nearest_station(pos);
        pos += dir;
        if (lock(pos)) { command(id, dir); }
    }

    void init() {
        int r = rows();
        int c = cols();
        for (int i = 0; i < r; ++i) {
            for (int j = 0; j < c; ++j) {
                Cell c = cell(i, j);
                if (c.type == Water) {
                    water_cells.insert(Pos(i, j));
                } else if (c.type == Station) {
                    station_cells.insert(Pos(i, j));
                }
            }
        }
        _nearest_water = _nearest_station =
            vector<vector<int>>(r, vector<int>(c, infinity));
        calculated_nearest_water = calculated_nearest_station = false;
        cell_locks = vector<vector<int>>(r, vector<int>(c, 0));
    }

public:
    void play() {
        if (round() == 0) { init(); }
        for (int car : cars(me())) {
            if (can_move(car)) { move_car(car); }
        }
        if (round() % 4 == me()) {
            for (int warrior : warriors(me())) { move_warrior(warrior); }
        }
    }

    static Player *factory() { return new PLAYER_NAME; }
};

RegisterPlayer(PLAYER_NAME);