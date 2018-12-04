#include "Player.hh"
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <utility>

#define PLAYER_NAME Linearithmic

using namespace std;

class PLAYER_NAME : public Player {
private:
    const int infinity = numeric_limits<int>::max();
    const int negative_infinity = numeric_limits<int>::min();
    const int number_of_players = 4;
    const unsigned int sentinels_per_city = 2;

    struct WarriorInfo {
        int city = -1;
        enum Role { Invader, TakingOver, Sentinel, Starving, Dehydrated };
        Role role = Invader;
    };
    map<int, WarriorInfo> warrior_info;

    /// Keeps a set of the city_cells for every city.
    vector<set<Pos>> city_cells;
    /// Contains a matrix of the distances to every city from every cell on the
    /// map.
    vector<vector<vector<int>>> nearest_city_;
    /// Contains a matrix of the distances to any water from every cell on the
    /// map.
    vector<vector<int>> nearest_water_;
    /// Contains a matrix of the **weighted** distances to any fuel station from
    /// every cell on the map.
    vector<vector<int>> nearest_station_;

    /// Stores the current owner of each city. Updated every moving turn.
    vector<int> city_owners;
    /// Marks if the owner changed, thus having to recalculate priorities.
    /// Updated every moving turn.
    vector<bool> city_owner_changed;
    /// Keeps track of the number of cities owned by each player.
    vector<int> city_count;
    /// Stores the number of warriors defined to attack an enemy city, or to
    /// assist a friendly one.
    vector<unsigned int> assigned_warriors;
    /// Stores the number of warriors defined to guard a friendly city.
    vector<unsigned int> assigned_sentinels;

    /** A vector containing the locked cells.
     * A lock is valid until the round marked by this vector (not included).
     * A lock may not be breached, since doing so may result in the death of
     * either of the units (or both!).
     */
    vector<vector<int>> cell_locks;

    enum UnitPlayer { Empty, Friend, Enemy, EnemyCar };

    inline UnitPlayer has_unit(int x, int y) {
        Cell c = cell(x, y);
        if (c.id == -1) { return Empty; }
        if (unit(c.id).player == me()) { return Friend; }
        if (unit(c.id).type == Car) { return EnemyCar; }
        return Enemy;
    }

    inline bool has_friend(int x, int y) { return has_unit(x, y) == Friend; }

    inline bool has_friend(Pos p) { return has_friend(p.i, p.j); }

    inline bool has_enemy(int x, int y, bool include_cars = true) {
        UnitPlayer up = has_unit(x, y);
        return up == Enemy or (include_cars and up == EnemyCar);
    }

    inline bool has_enemy(Pos p, bool include_cars = true) {
        return has_enemy(p.i, p.j, include_cars);
    }

    inline bool is_locked(int x, int y) {
        return has_friend(x, y) or cell_locks[x][y] > round();
    }

    inline bool is_locked(Pos p) { return is_locked(p.i, p.j); }

    inline bool lock(int x, int y) {
        if (pos_ok(x, y) and not is_locked(x, y)) {
            cell_locks[x][y] = round() + 1;
            return true;
        } else {
            return false;
        }
    }

    inline bool lock(Pos p) { return lock(p.i, p.j); }

    vector<Dir> random_dirs(bool include_none = false) {
        int upper_bound = include_none ? DirSize : (DirSize - 1);
        vector<Dir> ret(upper_bound);
        vector<int> rand_perm = random_permutation(upper_bound);
        for (unsigned int i = 0; i < ret.size(); ++i) {
            ret[i] = Dir(rand_perm[i]);
        }
        return ret;
    }

    Dir lock_or_escape(Pos p, Dir preferred_dir) {
        if (lock(p + preferred_dir)) { return preferred_dir; }
        if (lock(p)) { return None; }
        for (Dir d : random_dirs()) {
            if (lock(p + d)) { return d; }
        }
        return None;
    }

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
        _calculate_nearest(water_cells, nearest_water_, true, 0, 5);
    }

    inline void calculate_nearest_station(const set<Pos> &station_cells) {
        _calculate_nearest(station_cells, nearest_station_, false,
                           number_of_players - 1, 0);
    }

    inline void calculate_nearest_city(int city_id) {
        _calculate_nearest(city_cells[city_id], nearest_city_[city_id], true, 0,
                           5);
    }

    inline int _distance_to(Pos from, const vector<vector<int>> &distances) {
        return distances[from.i][from.j];
    }

    inline int distance_to_water(Pos from) {
        return _distance_to(from, nearest_water_);
    }

    inline int distance_to_station(Pos from) {
        return _distance_to(from, nearest_station_);
    }

    inline int distance_to_city(Pos from, int city_id) {
        return _distance_to(from, nearest_city_[city_id]);
    }

    inline int distance_to_nearest_city(Pos from) {
        int current_min = infinity;
        for (auto c : nearest_city_) {
            int dist_to_c = _distance_to(from, c);
            if (dist_to_c < current_min) { current_min = dist_to_c; }
        }
        return current_min;
    }

    Dir _nearest(Pos from, const vector<vector<int>> neighbourhood, const vector<vector<int>> &distances) {
        // TODO: Fix to use neighbourhood
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

    inline Dir nearest_water(Pos from, const vector<vector<int>> neighbourhood) {
        return _nearest(from, nearest_water_, neighbourhood);
    }

    inline Dir nearest_station(Pos from, const vector<vector<int>> neighbourhood) {
        return _nearest(from, nearest_station_, neighbourhood);
    }

    inline Dir goto_city(Pos from, int to_id, const vector<vector<int>> neighbourhood) {
        return _nearest(from, nearest_city_[to_id], neighbourhood);
    }

    Dir nearest_city(Pos from, const vector<vector<int>> &neighbourhood) {
        // TODO: Fix to use neighbourhood
        int current_min = infinity;
        int current_tgt = -1;
        for (unsigned int i = 0; i < nearest_city_.size(); ++i) {
            int dist_to_i = nearest_city_[i][from.i][from.j];
            if (dist_to_i < current_min) {
                current_min = dist_to_i;
                current_tgt = i;
            }
        }
        return (current_tgt == -1) ? None : goto_city(from, current_tgt);
    }

    void recalculate_city_owners() {
        vector<bool> changed(city_owner_changed.size(), false);
        for (unsigned int i = 0; i < city_cells.size(); ++i) {
            int old_owner = city_owners[i];
            int new_owner = cell(*(city_cells[i].begin())).owner;
            city_owners[i] = new_owner;
            changed[i] = old_owner != new_owner;
            if (changed[i] and old_owner != -1) {
                --city_count[old_owner];
                ++city_count[new_owner];
            }
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

    Dir dir_from_pos(Pos from, Pos to) {
        if (from.i < to.i) {
            if (from.j < to.j) {
                return Dir::BR;
            } else if (from.j > to.j) {
                return Dir::LB;
            } else { // from.j == to.j
                return Dir::Bottom;
            }
        } else if (from.i > to.i) {
            if (from.j < to.j) {
                return Dir::RT;
            } else if (from.j > to.j) {
                return Dir::TL;
            } else { // from.j == to.j
                return Dir::Top;
            }
        } else { // a.i == b.i
            if (from.j < to.j) {
                return Dir::Right;
            } else if (from.j > to.j) {
                return Dir::Left;
            } else { // from.j == to.j
                return Dir::None;
            }
        }
    }

    inline int score_warrior(const Unit &u) {
        if (u.type == Car) { return -10000; }
        return 10 + (min(u.food, u.water) + warriors(u.player).size() / 2 +
                     max(total_score(u.player) - total_score(me()), 0) * 4);
    }

    Dir find_nearest_enemy(Pos start_pos) {
        // Algorithm IDEA: Convert to Minimum(/Maximum) Spanning Tree and find
        // best path with DFS.
        queue<Pos> bfsq;
        bfsq.push(start_pos);
        set<Pos> visited;
        visited.insert(start_pos);
        Pos nearest;
        while (not bfsq.empty()) {
            Pos p = bfsq.front();
            bfsq.pop();
            Cell c = cell(p);
            if (c.type != Desert and c.type != Road) { continue; }
            if (c.id >= 0 and unit(c.id).player != me() and
                unit(c.id).type == Warrior) {
                nearest = p;
                break;
            }
            for (int d = 0; d < DirSize - 1; ++d) {
                Pos new_p = p + Dir(d);
                if (pos_ok(new_p) and visited.insert(new_p).second) {
                    bfsq.push(new_p);
                }
            }
        }
        return dir_from_pos(start_pos,
                            nearest); // A (bad) cheap approximation
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
        nearest_water_ = nearest_station_ =
            vector<vector<int>>(r, vector<int>(c, infinity));
        nearest_city_ = vector<vector<vector<int>>>(
            city_cells.size(),
            vector<vector<int>>(r, vector<int>(c, infinity)));
        calculate_nearest_water(water_cells);
        calculate_nearest_station(station_cells);
        for (unsigned int i = 0; i < nearest_city_.size(); ++i) {
            calculate_nearest_city(i);
        }
        cell_locks = vector<vector<int>>(r, vector<int>(c, 0));
        city_owners = vector<int>(city_cells.size());
        city_owner_changed = vector<bool>(city_owners.size(), -1);
        city_count = vector<int>(number_of_players, 0);
        assigned_warriors = assigned_sentinels =
            vector<unsigned int>(city_owners.size(), 0);
    }

    int score_city(Pos pos, int city_id) {
        int owner = city_owners[city_id];
        int dist = distance_to_city(pos, city_id);
        if (owner == me()) {
            if ((assigned_sentinels[city_id] + assigned_warriors[city_id]) <
                sentinels_per_city) {
                return 100 - dist;
            }
            return -dist * 10;
        }
        if (assigned_warriors[city_id] >= city_cells[city_id].size() / 3) {
            return (-dist * 3) / 2;
        }
        return 4 * city_count[owner] + total_score(owner) / 8 - dist;
    }

    Dir find_enemy_in_city(Pos start_pos, int limit) {
        queue<pair<Pos, int>> bfsq;
        bfsq.push({start_pos, 0});
        vector<Dir> directions = random_dirs();
        Pos preferred_pos = start_pos;
        Dir preferred_dir = None;
        while (not bfsq.empty()) {
            Pos current_pos = bfsq.front().first;
            int current_dist = bfsq.front().second;
            bfsq.pop();
            if (current_dist > limit or not pos_ok(current_pos) or
                is_locked(current_pos)) {
                continue;
            }
            Cell c = cell(current_pos);
            if (c.type != City) { continue; }
            if (has_enemy(current_pos, false)) {
                preferred_pos = current_pos;
                break;
            }
            for (Dir d : directions) {
                bfsq.push({current_pos + d, current_dist + 1});
            }
        }
        preferred_dir = dir_from_pos(start_pos, preferred_pos);
        if (preferred_dir == None) {
            for (Dir d : directions) {
                if (cell(start_pos + d).type == City) {
                    preferred_dir = d;
                    break;
                }
            }
        }
        return preferred_dir;
    }

    vector<vector<int>> look_around(Unit u) {
        constexpr int look_around_limit = 3;
        vector<vector<int>> neighbourhood(look_around_limit,
                                          vector<int>(look_around_limit, 0));
        Pos start_pos = u.pos;
        for (int i = 0; i < look_around_limit; ++i) {
            for (int j = 0; j < look_around_limit; ++i) {
                Pos curr_pos = start_pos + Pos(i, j);
                Cell c = cell(curr_pos);
                if (c.id == -1) { continue; }
                Unit u_in_cell = unit(c.id);
                if (u_in_cell.player == me()) { continue; }
                if (u_in_cell.type == Car) {
                    neighbourhood[i][j] -= 1000;
                    // TODO: Propagate to near cells
                }
                neighbourhood[i][j] +=
                    min(u_in_cell.food - u.food, u_in_cell.water - u.water);
            }
        }
        return neighbourhood;
    }

    inline bool check_water(const Unit &u, WarriorInfo &info,
                            unsigned int &assigned_to) {
        if (u.water < distance_to_water(u.pos) + 6) {
            info.role = WarriorInfo::Dehydrated;
            --assigned_to;
            return true;
        }
        return false;
    }

    inline bool check_food(const Unit &u, WarriorInfo &info,
                           unsigned int &assigned_to) {
        if (u.food < distance_to_nearest_city(u.pos) + 6) {
            info.role = WarriorInfo::Starving;
            --assigned_to;
            return true;
        }
        return false;
    }

    void move_warrior(int id, WarriorInfo &info) {
        Unit u = unit(id);
        Dir dir;
        switch (info.role) {
            case WarriorInfo::Starving: {
                if (u.food == warriors_health()) {
                    info.role = WarriorInfo::Invader;
                    info.city = -1;
                    return move_warrior(id, info);
                }
                vector<vector<int>> neighbourhood = move(look_around(u));
                dir = nearest_city(u.pos);
                break;
            }
            case WarriorInfo::Dehydrated: {
                if (u.water == warriors_health()) {
                    info.role = WarriorInfo::Invader;
                    info.city = -1;
                    return move_warrior(id, info);
                }
                vector<vector<int>> neighbourhood = look_around(u);
                dir = nearest_water(u.pos);
                break;
            }
            case WarriorInfo::Invader: {
                if (check_water(u, info, assigned_warriors[info.city]))
                    return move_warrior(id, info);
                if (info.city == -1 or city_owner_changed[info.city]) {
                    int current_max = negative_infinity;
                    int best_city = 0;
                    for (unsigned int city = 0; city < city_cells.size();
                         ++city) {
                        int score = score_city(u.pos, city);
                        if (score > current_max) {
                            current_max = score;
                            best_city = city;
                        }
                    }
                    info.city = best_city;
                    ++assigned_warriors[info.city];
                }
                vector<vector<int>> neighbourhood = move(look_around(u));
                if (distance_to_city(u.pos, info.city) == 0) {
                    info.role = WarriorInfo::TakingOver;
                    return move_warrior(id, info);
                }
                dir = goto_city(u.pos, info.city);
                break;
            }
            case WarriorInfo::TakingOver:
                if (check_water(u, info, assigned_warriors[info.city]))
                    return move_warrior(id, info);
                if (city_owners[info.city] == me()) {
                    info.role = WarriorInfo::Sentinel;
                    --assigned_warriors[info.city];
                    ++assigned_sentinels[info.city];
                }
                dir = find_enemy_in_city(u.pos, 5);
                break;
            case WarriorInfo::Sentinel:
                if (check_water(u, info, assigned_sentinels[info.city]))
                    return move_warrior(id, info);
                if (city_owners[info.city] != me()) {
                    info.role = WarriorInfo::TakingOver;
                    --assigned_sentinels[info.city];
                    ++assigned_warriors[info.city];
                }
                dir = find_enemy_in_city(u.pos, 2);
                break;
        }
        command(id, lock_or_escape(u.pos, dir));
    }

    inline void move_warrior(int id) { move_warrior(id, warrior_info[id]); }

    void move_car(int id) {
        Unit u = unit(id);
        Dir dir;
        if (distance_to_station(u.pos) >= u.food - 5) {
            dir = nearest_station(u.pos);
        } else {
            dir = find_nearest_enemy(u.pos);
        }
        if (lock(u.pos + dir)) {
            command(id, dir);
        } else if (not lock(u.pos)) {
            int d;
            for (d = 0; d < DirSize; ++d) {
                if (lock(u.pos + Dir(d))) {
                    command(id, Dir(d));
                    break;
                }
            }
            if (d == DirSize) { cerr << "Car " << id << " is doomed!" << endl; }
        }
    }

    void collect_dead_warriors() {
        const vector<int> &my_warriors = warriors(me());
        auto actual = my_warriors.begin();
        auto stored = warrior_info.begin();
        while (actual != my_warriors.end() and stored != warrior_info.end()) {
            if (*actual < stored->first) {
                ++actual;
            } else if (*actual > stored->first) {
                WarriorInfo &info = stored->second;
                if (info.city != -1 and
                    (info.role == WarriorInfo::Invader or
                     info.role == WarriorInfo::TakingOver)) {
                    --assigned_warriors[info.city];
                } else if (info.city != -1 and
                           info.role == WarriorInfo::Sentinel) {
                    --assigned_sentinels[info.city];
                }
                stored = warrior_info.erase(stored);
            } else {
                ++actual;
                ++stored;
            }
        }
        if (stored != warrior_info.end()) {
            warrior_info.erase(stored, warrior_info.end());
        }
    }

public:
    void play() {
        if (round() == 0) { init(); }
        for (int car : cars(me())) {
            if (can_move(car)) { move_car(car); }
        }
        if (round() % number_of_players == me()) {
            recalculate_city_owners();
            collect_dead_warriors();
            for (int warrior : warriors(me())) { move_warrior(warrior); }
        }
    }

    static Player *factory() { return new PLAYER_NAME; }
};

RegisterPlayer(PLAYER_NAME);