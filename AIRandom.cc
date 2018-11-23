#include "Player.hh"
#define PLAYER_NAME Random

using namespace std;

class PLAYER_NAME : public Player {
public:
    void play() {
        for (int id : cars(me())) {
            if (not can_move(id)) { continue; }
            command(id, Dir(random(0, DirSize - 1)));
        }
        if (round() % 4 != me()) { return; }
        for (int id : warriors(me())) {
            command(id, Dir(random(0, DirSize - 1)));
        }
    }

    /* Magic follows */
    static Player *factory() { return new PLAYER_NAME; }
};

/* More register magic */
RegisterPlayer(PLAYER_NAME);