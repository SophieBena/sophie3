//
//  BermanAlbite.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 1/19/11.
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "BermanProperties.h"

#define NS 2

@interface BermanAlbite : BermanProperties {
	@private
	double tOld, pOld, sOld[NS], invd2gds2[NS][NS];
	double gDis, hDis, sDis, cpDis, dcpdt, vDis, dvdt, dvdp, d2vdt2, d2vdtdp, d2vdp2;
}
@end

#undef NS
