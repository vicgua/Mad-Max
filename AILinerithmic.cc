#include "Player.hh"
#include <map>
#include <queue>

#define PLAYER_NAME Linearithmic

using namespace std;

/// Playing-class Linearithmic AI.
class PLAYER_NAME : public Player {
private:
    /** Possible Warrior modes. These indicate the current action done by the
     *  warrior.
     */
    enum WarriorMode {
        /// No task assigned. Random movements.
        Idle,
        /// Taking over a city.
        Invader,
        /// Protecting an owner city.
        Sentinel,
        /** Former Sentinel whose assigned city has been taking over.
         *  Will recover it or die.
         */
        Avenger,
        /** In need of supplies.
         *  Will go to the nearest city (if starving) or water tile (if
         *  dehydrated), avoiding enemies if possible (low chance of survival).
         */
        Survivor,
        /** Like Survivor, but unlikely to be able to resupply before dying.
         *  In a desperate move, will attack on sight any enemy warrior.
         */
        Doomed
    };

    /** Possible Car modes. These indicates the current action done by the car.
     */
    enum CarMode {
        /// No task assigned. Random movements.
        Wanderer,
        /// Out of fuel (if it needs fuel). Will go to the nearest Fuel Station.
        Refuel,
        /// Will chase enemy warriors and try to run over them.
        KillingSpree
    };

    /** Possible Car types. Unlike Modes, types are assigned on creation and
     *  remain unchanged until the unit is destroyed. */
    enum CarType {
        /// Will highly prefer roads over desert when moving.
        Patroller,
        /** Will go directly to the nearest enemy warrior, regardless of the
         *  tile type (unless city). Does not need to refuel. */
        OffRoad
    };

    struct Warrior {
        WarriorMode mode = Idle;
        Pos target = Pos(-1, -1);
    };

    struct Car {
        CarType type;
        CarMode mode = Wanderer;
    };

    /** A struct for assigning a score to a cell. */
    struct CellScore {
        Pos cell_pos;
        int score;
    };

    /** The AI is designed to operate in 4-times, so this will help in deciding
     * which tasks to do this round. */
    inline int time() const {
        return ((round() % 4) - me()) % 4;
        // Time 0: My turn.
    }

    static inline int dist_aprox(const Pos &from, const Pos &to) {
        return max(abs(from.i - to.i), abs(from.j - to.j));
    }

    static Dir dir_from_pos_dif(const Pos &from, const Pos &to) {
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

    map<int, Warrior> my_warriors;
    map<int, Car> my_cars;

    bool find_enemy(Pos start_pos, Pos &nearest_enemy, int max_depth,
                    bool can_city) const {
        queue<pair<Pos, int>> bfsq;
        bfsq.push({start_pos, max_depth});
        while (!bfsq.empty()) {
            auto curr = bfsq.front();
            bfsq.pop();
            if (not pos_ok(curr.first)) { continue; }
            Cell c = cell(curr.first);
            switch (c.type) {
                case Station:
                case Wall:
                case Water:
                    continue;
                case City:
                    if (not can_city) { continue; }
                    break;
                default:
                    break;
            }
            if (c.id != -1) {
                Unit u = unit(c.id);
                if (u.type == UnitType::Warrior and u.player != me()) {
                    nearest_enemy = curr.first;
                    return true;
                }
            }
            if (curr.second <= 0) { continue; } // Exceeded max_depth
            for (int d = 0; d < DirSize - 1; ++d) {
                bfsq.push({curr.first + Dir(d), curr.second - 1});
            }
        }
        return false;
    }

    void move_warrior(int id, Warrior &warrior_info) {
        // TODO
    }

    void move_car(int id, Car &car_info) {
        Unit unit_info = unit(id);
        if (car_info.type == Patroller) {
            // TODO
        } else if (car_info.type == OffRoad) {
            Pos dest_pos;
            Dir dest_dir;
            if (find_enemy(unit_info.pos, dest_pos, 3, false)) {
                dest_dir = dir_from_pos_dif(unit_info.pos, dest_pos);
            } else {
                dest_dir = Dir(random(0, DirSize - 1)); // None is not accepted
            }
            command(id, dest_dir);
        }
    }

    template <typename T>
    void recalculate(const vector<int> &from, map<int, T> &to) {
        queue<pair<typename map<int, T>::iterator, int>> captured;
        queue<typename map<int, T>::iterator> killed;
        auto actual = from.begin();
        auto stored = to.begin();
        while (actual != from.end() and stored != to.end()) {
            if (*actual < stored->first) {
                captured.push({stored, *actual});
                ++actual;
            } else if (*actual > stored->first) {
                killed.push(stored);
                ++stored;
            } else {
                ++actual;
                ++stored;
            }
        }
        while (not captured.empty()) {
            auto c = captured.front();
            captured.pop();
            to.insert(c.first, {c.second, T()});
        }
        while (not killed.empty()) {
            auto k = killed.front();
            killed.pop();
            to.erase(k);
        }
        while (actual != from.end()) {
            to.insert(to.end(), {*actual, T()});
            ++actual;
        }
        if (stored != to.end()) { to.erase(stored, to.end()); }
    }

    void recalculate_warriors() {
        recalculate<Warrior>(warriors(me()), my_warriors);
    }

    void recalculate_cars() { recalculate<Car>(cars(me()), my_cars); }

public:
    void play() {
        if (nb_players() != 4) {
            cerr << "Expected 4 players!" << endl;
            while (1)
                ;
        }
        int ctime = time();
        if (round() == 0 and ctime != 3) {
            recalculate_cars();
            recalculate_warriors();
        }

        switch (ctime) {
            case 0:
                for (auto &warrior : my_warriors) {
                    move_warrior(warrior.first, warrior.second);
                }
                break;
            case 1:
            case 2:
                break;
            case 3:
                recalculate_warriors();
                break;
            default:
                cerr << "Reached time " << time() << endl;
        }

        recalculate_cars();

        for (auto &car : my_cars) {
            if (can_move(car.first)) { move_car(car.first, car.second); }
        }
    }

    static Player *factory() { return new PLAYER_NAME; }
};

RegisterPlayer(PLAYER_NAME);