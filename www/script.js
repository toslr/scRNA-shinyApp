/**
 * @file script.js
 * @description Client-side JavaScript functionality for the shinyscRNA application
 * @author Tom Soulaire
 */

/**
 * Scrolls to a specific section with offset for the topbar
 * @param {string} sectionId - The ID of the section to scroll to
 */
function scrollToSection(sectionId) {
  const section = document.getElementById(sectionId);
  if (section) {
    const topbarHeight = document.getElementById('topbar').offsetHeight;
    const sectionRect = section.getBoundingClientRect();
    const scrollPosition = window.pageYOffset || document.documentElement.scrollTop;
    const targetPosition = scrollPosition + sectionRect.top - topbarHeight - 20; // 20px extra padding
    window.scrollTo({
      top: targetPosition,
      behavior: 'smooth'
    });
  }
}

// Initialize when document is ready
$(document).ready(function() {
  window.scrollToSection = scrollToSection;
  initializeCollapsibleSections();
  initializeCheckboxHandlers();
  setupDynamicElementObserver();
});

/**
 * Initialize collapsible sidebar sections
 */
function initializeCollapsibleSections() {
  // Set up click handlers for collapsible sections
  $('.collapsible-section .section-header').on('click', function() {
    var section = $(this).parent();
    section.toggleClass('collapsed');
    var sectionId = section.index();
    localStorage.setItem('section_' + sectionId + '_collapsed', section.hasClass('collapsed'));
  });
  
  // Initialize sections based on saved state (if any)
  $('.collapsible-section').each(function(index) {
    var defaultCollapsed = [2, 3, 4, 5, 6];
    var isCollapsed = defaultCollapsed.includes(index);
    if (isCollapsed) {
      $(this).addClass('collapsed');
    }
  });
}

/**
 * Initialize event handlers for sample, cluster, and condition checkboxes
 */
function initializeCheckboxHandlers() {
  // Sample management checkbox handler
  $(document).on('change', '.sample-management-checkbox', function() {
    var isChecked = $(this).prop('checked');
    var sampleId = $(this).attr('data-gsm');
    if (!sampleId) {
      var inputId = $(this).attr('id');
      sampleId = inputId.replace('sampleManagement-active_', '');
    }
    
    // Send the update to Shiny
    Shiny.setInputValue('sampleManagement-sample_checkbox_changed', {
      sample: sampleId,
      checked: isChecked
    });
  });
  
  // Cluster management checkbox handler
  $(document).on('change', '.cluster-management-checkbox', function() {
    var isChecked = $(this).prop('checked');
    var clusterId = $(this).attr('data-cluster');
    if (!clusterId) {
      var inputId = $(this).attr('id');
      clusterId = inputId.replace('clusterManagement-active_', '');
    }
    
    // Send the update to Shiny
    Shiny.setInputValue('clusterManagement-cluster_checkbox_changed', {
      cluster: clusterId,
      checked: isChecked
    });
  });
  
  // Condition management checkbox handler
  $(document).on('change', '.condition-management-checkbox', function() {
    var isChecked = $(this).prop('checked');
    var inputId = $(this).attr('id');
    
    // Delay the Shiny value update slightly to avoid race conditions
    setTimeout(function() {
      Shiny.setInputValue(inputId, isChecked);
    }, 10);
  });
}

/**
 * Setup a MutationObserver to monitor for dynamically added elements
 */
function setupDynamicElementObserver() {
  const observer = new MutationObserver(function(mutations) {
    mutations.forEach(function(mutation) {
      if (mutation.addedNodes && mutation.addedNodes.length > 0) {
        // Check if any condition checkboxes were added
        const conditionCheckboxes = $(mutation.addedNodes).find('input[id^="conditionManagement-active_"]');
        if (conditionCheckboxes.length > 0) {
          conditionCheckboxes.addClass('condition-management-checkbox');
        }
        
        // Check if any sample checkboxes were added
        const sampleCheckboxes = $(mutation.addedNodes).find('input[id^="sampleManagement-active_"]');
        if (sampleCheckboxes.length > 0) {
          sampleCheckboxes.addClass('sample-management-checkbox');
        }
        
        // Check if any cluster checkboxes were added
        const clusterCheckboxes = $(mutation.addedNodes).find('input[id^="clusterManagement-active_"]');
        if (clusterCheckboxes.length > 0) {
          clusterCheckboxes.addClass('cluster-management-checkbox');
        }
      }
    });
  });
  
  // Start observing the document body for changes
  observer.observe(document.body, { 
    childList: true, 
    subtree: true 
  });
}

// Custom message handlers
Shiny.addCustomMessageHandler("initializeSampleCheckboxes", function(message) {
  // This will be called when we need to initialize the sample checkboxes
  // Setup will be handled by the MutationObserver
});

Shiny.addCustomMessageHandler("initializeClusterCheckboxes", function(message) {
  // This will be called when we need to initialize the cluster checkboxes
  // Setup will be handled by the MutationObserver
});

// Handler for updating sample checkboxes
Shiny.addCustomMessageHandler("updateSampleCheckboxes", function(message) {
  const selectedSamples = message.samples;
  
  // Find all sample checkboxes
  const checkboxes = document.querySelectorAll('input.sample-management-checkbox');
  
  // Update each checkbox
  checkboxes.forEach(function(checkbox) {
    let sampleId = checkbox.getAttribute('data-gsm');
    if (!sampleId) {
      const inputId = checkbox.getAttribute('id');
      if (inputId) {
        sampleId = inputId.replace('sampleManagement-active_', '');
      }
    }
    if (sampleId) {
      checkbox.checked = selectedSamples.includes(sampleId);
    }
  });
  
  // Also update the "select all" checkbox
  const selectAllCheckbox = document.getElementById('sampleManagement-selectAllSamples');
  if (selectAllCheckbox) {
    selectAllCheckbox.checked = checkboxes.length > 0 && 
      Array.from(checkboxes).every(cb => cb.checked);
  }
});

// Handler for updating cluster checkboxes
Shiny.addCustomMessageHandler("updateClusterCheckboxes", function(message) {
  const selectedClusters = message.clusters;
  
  // Find all cluster checkboxes
  const checkboxes = document.querySelectorAll('input.cluster-management-checkbox');
  
  // Update each checkbox
  checkboxes.forEach(function(checkbox) {
    let clusterId = checkbox.getAttribute('data-cluster');
    if (!clusterId) {
      const inputId = checkbox.getAttribute('id');
      if (inputId) {
        clusterId = inputId.replace('clusterManagement-active_', '');
      }
    }
    if (clusterId) {
      checkbox.checked = selectedClusters.includes(clusterId);
    }
  });
  
  // Also update the "select all" checkbox
  const selectAllCheckbox = document.getElementById('clusterManagement-selectAllClusters');
  if (selectAllCheckbox) {
    selectAllCheckbox.checked = checkboxes.length > 0 && 
      Array.from(checkboxes).every(cb => cb.checked);
  }
});