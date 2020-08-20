//
// Created by carlos on 06/03/19.
//

#ifndef MRP_ARC_H
#define MRP_ARC_H


class Arc {

private:
    int o, d;
    int delay, jitter, bandwidth, estimateLinkDuration;

public:
    Arc(int o, int d, int delay, int jitter, int bandwidth, int estimateLinkDuration);

    int getO() const;

    int getD() const;

    int getDelay() const;

    int getJitter() const;

    int getBandwidth() const;

    int getEstimateLinkDuration() const;
};


#endif //MRP_ARC_H
