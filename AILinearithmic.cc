#include "Player.hh"
#include <limits>
#include <queue>
#include <set>
#include <utility>
#include <functional>

#define PLAYER_NAME Linearithmic

using namespace std;

class PLAYER_NAME : public Player {
private:
    const int infinity = numeric_limits<int>::max();

    vector<set<Pos>> city_cells;
    vector<vector<int>> _nearest_water;
    vector<vector<int>> _nearest_station;

    vector<int> city_owners;
    vector<bool> city_owner_changed;

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

    void _calculate_nearest(const set<Pos> &cells,
                            vector<vector<int>> &dist_matrix,
                            const bool can_pass_city, const int desert_penalty,
                            const int road_penalty) {
        queue<Pos> bfsq;
        for (auto pos : cells) {
            dist_matrix[pos.i][pos.j] = 0;
            bfsq.push(pos);
        }
        while (not bfsq.empty()) {
            Pos curr = bfsq.front();
            bfsq.pop();
            int base_dist = dist_matrix[curr.i][curr.j];
            for (int d = 0; d < DirSize; ++d) {
                int new_dist = base_dist + 1;
                Pos next_pos = curr + Dir(d);
                if (not pos_ok(next_pos)) { continue; }
                switch (cell(next_pos).type) {
                    case Station:
                    case Wall:
                    case Water:
                        continue; // Impassable
                    case City:
                        if (not can_pass_city) { continue; }
                        break;
                    case Desert:
                        new_dist += desert_penalty;
                        break;
                    case Road:
                        new_dist += road_penalty;
                        break;
                    default:
                        // TODO: Assert unreachable
                        break;
                }
                if (dist_matrix[next_pos.i][next_pos.j] > new_dist) {
                    dist_matrix[next_pos.i][next_pos.j] = new_dist;
                    bfsq.push(next_pos);
                }
            }
        }
    }

    inline void calculate_nearest_water(const set<Pos> &water_cells) {
        _calculate_nearest(water_cells, _nearest_water, true, 0, 1);
    }

    inline void calculate_nearest_station(const set<Pos> &station_cells) {
        _calculate_nearest(station_cells, _nearest_station, false, 3, 0);
    }

    Dir _nearest(Pos from, const vector<vector<int>> &distances) {
        int current_min = distances[from.i][from.j];
        Dir current_dir = None;
        vector<int> directions(random_permutation(DirSize - 1));
        for (int d : directions) {
            Pos next_pos = from + Dir(d);
            if (not pos_ok(next_pos)) { continue; }
            int next_pos_dist = distances[next_pos.i][next_pos.j];
            if (next_pos_dist < current_min) {
                current_min = next_pos_dist;
                current_dir = Dir(d);
            }
        }
        return current_dir;
    }

    inline Dir nearest_water(Pos from) {
        return _nearest(from, _nearest_water);
    }

    inline Dir nearest_station(Pos from) {
        return _nearest(from, _nearest_station);
    }

    void recalculate_city_owners() {
        vector<bool> changed(city_owner_changed.size(), false);
        for (unsigned int i = 0; i < city_cells.size(); ++i) {
            int old_owner = city_owners[i];
            int new_owner = cell(*(city_cells[i].begin())).owner;
            city_owners[i] = new_owner;
            changed[i] = old_owner != new_owner;
        }
        city_owner_changed = move(changed);
    }

    void add_city_cell(vector<vector<int>> &city_map, Pos pos) {
        for (int d = 0; d < DirSize - 1; ++d) {
            Pos possible_pos = pos + Dir(d);
            if (not pos_ok(possible_pos)) { continue; }
            int possible_city = city_map[possible_pos.i][possible_pos.j];
            if (possible_city >= 0) {
                city_map[pos.i][pos.j] = possible_city;
                city_cells[possible_city].insert(pos);
                return;
            }
        }
        int new_city_id = city_cells.size();
        city_map[pos.i][pos.j] = new_city_id;
        set<Pos> new_city{pos};
        city_cells.push_back(move(new_city));
    }

    void init() {
        int r = rows();
        int c = cols();
        vector<vector<int>> city_map(r, vector<int>(c, -1));
        set<Pos> water_cells;
        set<Pos> station_cells;
        for (int i = 0; i < r; ++i) {
            for (int j = 0; j < c; ++j) {
                Cell c = cell(i, j);
                switch (c.type) {
                    case City:
                        add_city_cell(city_map, Pos(i, j));
                        break;
                    case Station:
                        station_cells.insert(Pos(i, j));
                        break;
                    case Water:
                        water_cells.insert(Pos(i, j));
                        break;
                    default:
                        break;
                }
            }
        }
        _nearest_water = _nearest_station =
            vector<vector<int>>(r, vector<int>(c, infinity));
        calculate_nearest_water(water_cells);
        calculate_nearest_station(station_cells);
        cell_locks = vector<vector<int>>(r, vector<int>(c, 0));
        city_owners = vector<int>(city_cells.size());
        city_owner_changed = vector<bool>(city_owners.size());
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


public:
    void play() {
        if (round() == 0) { init(); }
        for (int car : cars(me())) {
            if (can_move(car)) { move_car(car); }
        }
        if (round() % 4 == me()) {
            recalculate_city_owners();
            for (int warrior : warriors(me())) { move_warrior(warrior); }
        }
    }

    static Player *factory() { return new PLAYER_NAME; }
};

RegisterPlayer(PLAYER_NAME);