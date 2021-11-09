class BoundException(Exception):
    def __init__(self, participants, message="clipped with each other"):
        self.participants = participants
        self.message = message
        string = 'Participants: {} {}'.format(participants, self.message)
        super(BoundException, self).__init__(string)


class UnboundException(Exception):
    def __init__(self, particle, message="has escaped all bounds."):
        self.particle = particle
        self.message = message
        string = 'Particle: {} {}'.format(particle, self.message)
        super(UnboundException, self).__init__(string)
