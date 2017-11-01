class GlobalData:
    def __init__(self):
        self.spectral_region = WvlSpec()
        self.pupil = PupilSpec()
        self.field_of_view = FieldSpec()
        self.specs = SystemSpec()


class WvlSpec:
    def __init__(self):
        self.spectrum = []
        self.reference_wvl = 0
        self.coating_wvl = 550.0

    def add(self, wl, wt):
        self.spectrum.append([wl, wt])
        self.spectrum.sort(key=lambda w: w[0], reverse=True)


class PupilSpec:
    types = ('EPD', 'NA', 'NAO', 'FNO')

    def __init__(self):
        self.type = 'EPD'
        self.value = 1.0


class FieldSpec:
    types = ('OBJ_ANG', 'OBJ_HT', 'IMG_HT')

    def __init__(self):
        self.fields = []
        self.type = 'OBJ_ANG'
        self.wide_angle = False


class Field:
    def __init__(self):
        self.x = 0.0
        self.y = 0.0
        self.vux = 0.0
        self.vuy = 0.0
        self.vlx = 0.0
        self.vly = 0.0
        self.wt = 1.0


class SystemSpec:
    dims = ('M', 'C', 'I')

    def __init__(self):
        self.radius_mode = False
        self.title = ''
        self.initials = ''
        self.dimensions = 'M'
        self.aperture_override = ''
        self.temperature = 20.0
        self.pressure = 760.0
