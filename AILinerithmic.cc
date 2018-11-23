#include "Player.hh"

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
        WarriorMode mode = Wanderer;
        Pos target = Pos(-1, -1);
    };

    struct Car {
        CarType type;
        CarMode mode;
    };

public:
    void play() {}

    static Player *factory() { return new PLAYER_NAME; }
};

RegisterPlayer(PLAYER_NAME);