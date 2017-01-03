import json
import pvl
import spiceypy as spice

def find_in_dict(obj, key):
    """
    Recursively find an entry in a dictionary

    Parameters
    ----------
    obj : dict
          The dictionary to search
    key : str
          The key to find in the dictionary

    Returns
    -------
    item : obj
           The value from the dictionary
    """
    if key in obj:
        return obj[key]
    for k, v in obj.items():
        if isinstance(v,dict):
            item = find_in_dict(v, key)
            if item is not None:
                return item


def main(kernelid=236820):
    spice.furnsh("../tests/data/msgr_mdis_v160.ti")
    isd  = {}

    # Load information from the IK kernel
    isd['focal_length'] = spice.gdpool('INS-{}_FOCAL_LENGTH'.format(kernelid), 0, 1).tolist()[0]
    isd['focal_length_epsilon'] = spice.gdpool('INS-{}_FL_UNCERTAINTY'.format(kernelid), 0, 1).tolist()[0]
    isd['nlines'] = spice.gipool('INS-{}_PIXEL_LINES'.format(kernelid), 0, 1).tolist()[0]
    isd['nsamples'] = spice.gipool('INS-{}_PIXEL_SAMPLES'.format(kernelid), 0, 1).tolist()[0]
    isd['original_half_lines'] = isd['nlines'] / 2.0
    isd['orginal_half_samples'] = isd['nsamples'] / 2.0
    isd['pixel_pitch'] = spice.gdpool('INS-{}_PIXEL_PITCH'.format(kernelid), 0, 1).tolist()[0]
    isd['ccd_center'] = spice.gdpool('INS-{}_CCD_CENTER'.format(kernelid), 0, 1).tolist()[0]
    isd['ifov'] = spice.gdpool('INS-{}_IFOV'.format(kernelid), 0, 1).tolist()[0]
    isd['boresight'] = spice.gdpool('INS-{}_BORESIGHT'.format(kernelid), 0, 3).tolist()
    isd['transx'] = spice.gdpool('INS-{}_TRANSX'.format(kernelid), 0, 3).tolist()
    isd['transy'] = spice.gdpool('INS-{}_TRANSY'.format(kernelid), 0, 3).tolist()
    isd['itrans_sample'] = spice.gdpool('INS-{}_ITRANSS'.format(kernelid), 0, 3).tolist()[0]
    isd['itrans_line'] = spice.gdpool('INS-{}_ITRANSL'.format(kernelid), 0, 3).tolist()[0]
    isd['odt_x'] = spice.gdpool('INS-{}_OD_T_X'.format(kernelid), 0, 9).tolist()
    isd['odt_y'] = spice.gdpool('INS-{}_OD_T_Y'.format(kernelid), 0, 9).tolist()
    isd['starting_detector_sample'] = spice.gdpool('INS-{}_FPUBIN_START_SAMPLE'.format(kernelid), 0, 1).tolist()[0]
    isd['starting_detector_line'] = spice.gdpool('INS-{}_FPUBIN_START_LINE'.format(kernelid), 0, 1).tolist()[0]

    # Load the ISIS Cube header
    header = pvl.load('../tests/data/CN0108840044M_IF_5_NAC_spiced.cub')

    isd['instrument_id'] = find_in_dict(header, 'InstrumentId')
    isd['spacecraft_name'] = find_in_dict(header, 'SpacecraftName')

    #Time
    # Load LeapSecond Kernel
    spice.furnsh('../tests/data/naif0011.tls.txt')
    start_ephemeris_time = spice.str2et(find_in_dict(header, 'StartTime').isoformat())
    stop_ephemeris_time = spice.str2et(find_in_dict(header, 'StopTime').isoformat())

    # Total hack - short on time - set the et to be the start time - this is WRONG
    isd['ephemeris_time'] = start_ephemeris_time

    # OPK and Position - the logic to get these is in Anne's SocetCode and the
    # called ISIS functions - starts on line 418 of the cpp
    isd['x_sensor_origin'] = None
    isd['y_sensor_origin'] = None
    isd['z_sensor_origin'] = None
    isd['omega'] = None
    isd['phi'] = None
    isd['kappa'] = None

    # ISD Search Information - totally fabricated - where do we get these?
    isd['min_elevation'] = -1.0
    isd['max_elevation'] = 1.0

    # Write it out
    with open('hardcoded.isd', 'w') as f:
        json.dump(isd, f, sort_keys=True, indent=4)

if __name__ == '__main__':
    main()
