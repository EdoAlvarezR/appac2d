import os, shutil

def create_path(save_path, prompt):
    
    # Create/flush save path
    if save_path != None:
    
        create = True
        
        if os.path.exists(save_path):
    
            response = None if prompt else 'y'
            while response not in ['y', 'n']:
                response = input(f"Folder '{save_path}' already exists. Remove? (y/n) ").strip().lower()
    
            if response == 'y':
                shutil.rmtree(save_path)
            else:
                create = False
    
        if create:
            os.makedirs(save_path)


def PR_to_CT(atmosphere, PR=1.0, altitude_ft=0.0, Mach=0.2):

    altitude            = 0.3048 * altitude_ft           # (m)
    
    Patm = atmosphere.pressure(altitude)
    Uinf = atmosphere.speed_of_sound(altitude) * Mach
    rho = atmosphere.density(altitude)
        
    P0 = Patm + 0.5 * rho * Uinf**2
    P1 = PR * P0

    deltaP = P1 - P0

    CT = deltaP / (rho*Uinf**2)
        
    return CT