from  chimera.extension import EMO, manager

# -----------------------------------------------------------------------------
#
class MapQ_Dialog_EMO ( EMO ):

  def name(self):
    return 'MapQ'
  def description(self):
    return self.categoryDescriptions()['Volume Data']
  def categories(self):
    return self.categoryDescriptions().keys()
  def categoryDescriptions(self):
    # since we want to use specialized descriptions for certain categories...
    return {
      'Volume Data': 'Evaluate map & model',
    }
  def icon(self):
    return self.path('mapq.png')
  def activate(self):
    # self.module('volumedialog').show_volume_dialog()
    d = self.module('mapq').show_dialog()
    return None

# -----------------------------------------------------------------------------
# Register dialogs and menu entry.
#
manager.registerExtension ( MapQ_Dialog_EMO ( __file__ ) )
